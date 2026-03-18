#!/usr/bin/env python
"""
FDTD (Meep) validation script for miepy's C++ T-matrix solver.

Runs Meep scattering simulations and compares cross-sections against miepy.
Results are cached to disk to avoid re-running expensive FDTD simulations.

Usage:
    pixi run python tests/fdtd_validation.py --list
    pixi run python tests/fdtd_validation.py --run sphere_dielectric
    pixi run python tests/fdtd_validation.py --run all
    pixi run python tests/fdtd_validation.py --run sphere_dielectric --plot
    pixi run python tests/fdtd_validation.py --run sphere_dielectric --force
    pixi run python tests/fdtd_validation.py --run all --resolution high
    pixi run python tests/fdtd_validation.py --summary
"""

import argparse
import sys
from pathlib import Path

import numpy as np

nm = 1e-9
um = 1e-6

CACHE_DIR = Path(__file__).parent / "fdtd_cache"

RESOLUTION_PRESETS = {
    "fast": 1 / (12 * nm),
    "standard": 1 / (8 * nm),
    "high": 1 / (5 * nm),
    "veryhigh": 1 / (4 * nm),
    "extreme": 1 / (2 * nm),
}

# ---------------------------------------------------------------------------
# Inline helpers (replacing meep_ext)
# ---------------------------------------------------------------------------


def freq_data(fmin, fmax):
    """Return (fcen, df) for a GaussianSource spanning [fmin, fmax]."""
    fcen = 0.5 * (fmin + fmax)
    df = fmax - fmin
    return fcen, df


def add_flux_box(sim, fcen, df, nfreq, center, size):
    """Add 6 FluxRegion faces for net outward flux through a box."""
    import meep

    cx, cy, cz = center
    sx, sy, sz = size
    V = meep.Vector3

    regions = [
        meep.FluxRegion(center=V(cx + sx / 2, cy, cz), size=V(0, sy, sz), weight=+1),
        meep.FluxRegion(center=V(cx - sx / 2, cy, cz), size=V(0, sy, sz), weight=-1),
        meep.FluxRegion(center=V(cx, cy + sy / 2, cz), size=V(sx, 0, sz), weight=+1),
        meep.FluxRegion(center=V(cx, cy - sy / 2, cz), size=V(sx, 0, sz), weight=-1),
        meep.FluxRegion(center=V(cx, cy, cz + sz / 2), size=V(sx, sy, 0), weight=+1),
        meep.FluxRegion(center=V(cx, cy, cz - sz / 2), size=V(sx, sy, 0), weight=-1),
    ]

    return sim.add_flux(fcen, df, nfreq, *regions)


def add_flux_plane(sim, fcen, df, nfreq, center, size):
    """Add a single planar flux monitor for incident power normalization."""
    import meep

    cx, cy, cz = center
    sx, sy, sz = size
    region = meep.FluxRegion(
        center=meep.Vector3(cx, cy, cz),
        size=meep.Vector3(sx, sy, sz),
    )
    return sim.add_flux(fcen, df, nfreq, region)


def make_plane_wave_sources(src_time, cell_size, pml_thickness):
    """Create x-polarized plane wave sources at the bottom of the cell."""
    import meep

    src_z = -0.5 * cell_size.z + pml_thickness
    return [
        meep.Source(
            src_time,
            component=meep.Ex,
            center=meep.Vector3(0, 0, src_z),
            size=meep.Vector3(cell_size.x, cell_size.y, 0),
        )
    ]


# ---------------------------------------------------------------------------
# Meep scattering simulation
# ---------------------------------------------------------------------------


class MeepScatteringSim:
    """Encapsulates the two-simulation FDTD scattering pattern.

    Args:
        geometry: Meep geometry objects (in Meep coordinate units).
        box_size: Particle bounding box [x, y, z] in meters.
        resolution: Grid resolution in pixels per meter.
        wl_min, wl_max: Wavelength range in meters.
        nfreq: Number of frequency points for DFT monitors.
        medium_index: Background refractive index.
        length_scale: Meep length unit in meters. Default 1.0 (SI meters).
            Set to 1e-6 for um-based coordinates (needed for Meep Prism geometry,
            which has numerical issues at SI-meter coordinate scales).
            When using length_scale != 1.0, geometry must be pre-scaled to Meep units.
    """

    def __init__(self, geometry, box_size, resolution, wl_min, wl_max,
                 nfreq=40, medium_index=1.0, length_scale=1.0):
        import meep

        self.geometry = geometry
        self.nfreq = nfreq
        self.medium = meep.Medium(index=medium_index)
        self.length_scale = length_scale
        ls = length_scale

        # Convert resolution: pixels/m -> pixels/Meep_unit
        self.resolution = resolution * ls

        # Gaps in Meep units
        particle_monitor_gap = 50 * nm / ls
        pml_monitor_gap = 50 * nm / ls
        self.pml_thickness = 15 / self.resolution

        # Box size in Meep units
        box_meep = [b / ls for b in box_size]

        # Monitor box around particle (Meep units)
        self.monitor_size = [b + 2 * particle_monitor_gap for b in box_meep]

        # Cell size (Meep units)
        self.cell = meep.Vector3(
            *[b + 2 * particle_monitor_gap + 2 * pml_monitor_gap + 2 * self.pml_thickness for b in box_meep]
        )

        self.pml = [meep.PML(self.pml_thickness)]

        # Flux monitor frequency range (cycles per Meep unit)
        self.fcen, self.df = freq_data(ls / wl_max, ls / wl_min)

        # Source: match measurement bandwidth for maximum signal-to-noise.
        self.src_time = meep.GaussianSource(frequency=self.fcen, fwidth=self.df)

        # Area for cross-section normalization (Meep units^2)
        self.area_meep = box_meep[0] * box_meep[1]
        self.box_meep = box_meep

    def _make_sources(self):
        return make_plane_wave_sources(self.src_time, self.cell, self.pml_thickness)

    def _run_norm(self):
        import meep

        sim = meep.Simulation(
            cell_size=self.cell,
            boundary_layers=self.pml,
            geometry=[],
            default_material=self.medium,
            resolution=self.resolution,
            sources=self._make_sources(),
        )

        flux_box = add_flux_box(sim, self.fcen, self.df, self.nfreq, [0, 0, 0], self.monitor_size)
        flux_plane = add_flux_plane(
            sim, self.fcen, self.df, self.nfreq, [0, 0, 0], [self.box_meep[0], self.box_meep[1], 0]
        )

        # decay dt in Meep units (0.5 um equivalent)
        decay_dt = 0.5 * um / self.length_scale
        sim.run(
            until_after_sources=meep.stop_when_fields_decayed(
                decay_dt,
                meep.Ex,
                pt=meep.Vector3(0, 0, self.monitor_size[2] / 2),
                decay_by=1e-3,
            )
        )

        norm_flux_data = sim.get_flux_data(flux_box)

        result = {
            "frequency": np.array(meep.get_flux_freqs(flux_plane)),
            "incident": np.array(meep.get_fluxes(flux_plane)),
        }

        scat_result = self._run_scat(norm_flux_data)

        return result, scat_result

    def _run_scat(self, norm_flux_data):
        import meep

        sim = meep.Simulation(
            cell_size=self.cell,
            boundary_layers=self.pml,
            geometry=self.geometry,
            default_material=self.medium,
            resolution=self.resolution,
            sources=self._make_sources(),
        )

        flux_box_absorb = add_flux_box(sim, self.fcen, self.df, self.nfreq, [0, 0, 0], self.monitor_size)
        flux_box_scat = add_flux_box(sim, self.fcen, self.df, self.nfreq, [0, 0, 0], self.monitor_size)
        sim.load_minus_flux_data(flux_box_scat, norm_flux_data)

        decay_dt = 0.5 * um / self.length_scale
        sim.run(
            until_after_sources=meep.stop_when_fields_decayed(
                decay_dt,
                meep.Ex,
                pt=meep.Vector3(0, 0, self.monitor_size[2] / 2),
                decay_by=1e-5,
            )
        )

        return {
            "scattering": np.array(meep.get_fluxes(flux_box_scat)),
            "absorption": -np.array(meep.get_fluxes(flux_box_absorb)),
            "frequency": np.array(meep.get_flux_freqs(flux_box_scat)),
        }

    def run(self):
        """Run norm + scat simulations and return cross-section data."""
        norm_data, scat_data = self._run_norm()

        frequency = norm_data["frequency"]
        incident = norm_data["incident"]

        # Cross-sections in Meep units^2, convert to m^2
        ls2 = self.length_scale**2
        cs_scat = scat_data["scattering"] / incident * self.area_meep * ls2
        cs_abs = scat_data["absorption"] / incident * self.area_meep * ls2
        cs_ext = cs_scat + cs_abs

        # Frequencies in cycles/Meep_unit, convert to cycles/m
        frequency_si = frequency / self.length_scale

        return {
            "frequency": frequency_si,
            "cross_section_scat": cs_scat,
            "cross_section_abs": cs_abs,
            "cross_section_ext": cs_ext,
        }


class MeepScatteringSimGold:
    """Meep scattering simulation using um-based coordinates for gold (meep.materials.Au)."""

    def __init__(self, geometry, box_size_um, resolution_um, wl_min_um, wl_max_um, nfreq=40):
        import meep

        self.geometry = geometry
        self.resolution = resolution_um
        self.nfreq = nfreq
        self.medium = meep.Medium(index=1)

        particle_monitor_gap = 0.05  # 50 nm in um
        pml_monitor_gap = 0.05
        pml_thickness = 15 / resolution_um

        self.monitor_size = [b + 2 * particle_monitor_gap for b in box_size_um]

        self.cell = meep.Vector3(
            *[b + 2 * particle_monitor_gap + 2 * pml_monitor_gap + 2 * pml_thickness for b in box_size_um]
        )

        self.pml = [meep.PML(pml_thickness)]
        self.fcen, self.df = freq_data(1 / wl_max_um, 1 / wl_min_um)
        self.src_time = meep.GaussianSource(frequency=self.fcen, fwidth=self.df)
        self.area_um2 = box_size_um[0] * box_size_um[1]
        self.box_size_um = box_size_um
        self.pml_thickness = pml_thickness

    def _make_sources(self):
        return make_plane_wave_sources(self.src_time, self.cell, self.pml_thickness)

    def _run_norm(self):
        import meep

        sim = meep.Simulation(
            cell_size=self.cell,
            boundary_layers=self.pml,
            geometry=[],
            default_material=self.medium,
            resolution=self.resolution,
            sources=self._make_sources(),
        )

        flux_box = add_flux_box(sim, self.fcen, self.df, self.nfreq, [0, 0, 0], self.monitor_size)
        flux_plane = add_flux_plane(
            sim, self.fcen, self.df, self.nfreq, [0, 0, 0], [self.box_size_um[0], self.box_size_um[1], 0]
        )

        sim.run(
            until_after_sources=meep.stop_when_fields_decayed(
                0.5,
                meep.Ex,
                pt=meep.Vector3(0, 0, self.monitor_size[2] / 2),
                decay_by=1e-3,
            )
        )

        norm_flux_data = sim.get_flux_data(flux_box)

        norm_data = {
            "frequency": np.array(meep.get_flux_freqs(flux_plane)),
            "incident": np.array(meep.get_fluxes(flux_plane)),
        }

        scat_data = self._run_scat(norm_flux_data)

        return norm_data, scat_data

    def _run_scat(self, norm_flux_data):
        import meep

        sim = meep.Simulation(
            cell_size=self.cell,
            boundary_layers=self.pml,
            geometry=self.geometry,
            default_material=self.medium,
            resolution=self.resolution,
            sources=self._make_sources(),
        )

        flux_box_absorb = add_flux_box(sim, self.fcen, self.df, self.nfreq, [0, 0, 0], self.monitor_size)
        flux_box_scat = add_flux_box(sim, self.fcen, self.df, self.nfreq, [0, 0, 0], self.monitor_size)
        sim.load_minus_flux_data(flux_box_scat, norm_flux_data)

        sim.run(
            until_after_sources=meep.stop_when_fields_decayed(
                0.5,
                meep.Ex,
                pt=meep.Vector3(0, 0, self.monitor_size[2] / 2),
                decay_by=1e-5,
            )
        )

        return {
            "scattering": np.array(meep.get_fluxes(flux_box_scat)),
            "absorption": -np.array(meep.get_fluxes(flux_box_absorb)),
            "frequency": np.array(meep.get_flux_freqs(flux_box_scat)),
        }

    def run(self):
        """Run norm + scat simulations, return cross-sections in m^2."""
        norm_data, scat_data = self._run_norm()

        frequency = norm_data["frequency"]
        incident = norm_data["incident"]

        # Cross-sections in um^2, convert to m^2
        cs_scat = scat_data["scattering"] / incident * self.area_um2 * um**2
        cs_abs = scat_data["absorption"] / incident * self.area_um2 * um**2
        cs_ext = cs_scat + cs_abs

        # Frequencies are in cycles/um, convert to cycles/m for consistency
        frequency_si = frequency / um

        return {
            "frequency": frequency_si,
            "cross_section_scat": cs_scat,
            "cross_section_abs": cs_abs,
            "cross_section_ext": cs_ext,
        }


# ---------------------------------------------------------------------------
# miepy cross-section computation
# ---------------------------------------------------------------------------


def compute_miepy_cross_sections(particle, wavelengths, lmax, medium=None):
    """Compute scattering, absorption, extinction cross-sections via miepy."""
    import miepy

    source = miepy.sources.plane_wave([1, 0])
    C_scat = np.zeros(len(wavelengths))
    C_abs = np.zeros(len(wavelengths))
    C_ext = np.zeros(len(wavelengths))

    for i, wl in enumerate(wavelengths):
        sol = miepy.cluster(
            particles=[particle],
            source=source,
            wavelength=wl,
            lmax=lmax,
            medium=medium,
        )
        C_scat[i], C_abs[i], C_ext[i] = sol.cross_sections()

    return C_scat, C_abs, C_ext


# ---------------------------------------------------------------------------
# Disk caching
# ---------------------------------------------------------------------------


def cache_path(test_name, resolution):
    """Get cache file path for a test case."""
    res_nm = round(1 / (resolution * nm))
    return CACHE_DIR / f"{test_name}_res{res_nm}nm.npz"


def load_cache(test_name, resolution):
    """Load cached Meep results if available."""
    path = cache_path(test_name, resolution)
    if path.exists():
        data = np.load(path)
        return {k: data[k] for k in data.files}
    return None


def save_cache(test_name, resolution, data):
    """Save Meep results to disk."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = cache_path(test_name, resolution)
    np.savez(path, **data)
    print(f"  Cached to {path}")


def load_or_run(test_name, resolution, run_func, force=False):
    """Load from cache or run the simulation."""
    if not force:
        cached = load_cache(test_name, resolution)
        if cached is not None:
            print(f"  Loaded from cache")
            return cached

    print(f"  Running Meep simulation...")
    data = run_func(resolution)
    save_cache(test_name, resolution, data)
    return data


# ---------------------------------------------------------------------------
# Comparison metrics
# ---------------------------------------------------------------------------


def compute_errors(meep_data, miepy_data, meep_freqs):
    """Compute relative errors between Meep and miepy cross-sections.

    miepy_data is at 100 dense wavelength points; interpolate to Meep frequencies.
    """
    miepy_wavelengths = miepy_data["wavelengths"]
    miepy_freqs = 1.0 / miepy_wavelengths  # cycles/m

    # Use extinction as reference scale for all error denominators
    meep_ext = meep_data["cross_section_ext"]
    ref_scale = np.max(np.abs(meep_ext))

    results = {}
    for cs_type in ["scat", "abs", "ext"]:
        meep_cs = meep_data[f"cross_section_{cs_type}"]
        miepy_cs_key = f"cross_section_{cs_type}"

        # Interpolate miepy to Meep frequency points
        miepy_cs = np.interp(meep_freqs, miepy_freqs[::-1], miepy_data[miepy_cs_key][::-1])

        # Use max of own signal or small fraction of extinction as denominator
        # This avoids division by zero for absorption of dielectrics
        denom = np.abs(meep_cs)
        threshold = ref_scale * 1e-3
        denom = np.where(denom > threshold, denom, threshold)
        rel_err = np.abs(meep_cs - miepy_cs) / denom

        meep_norm = np.sum(meep_cs**2)
        if meep_norm > (ref_scale * 1e-3) ** 2 * len(meep_cs):
            rel_L2 = np.sqrt(np.sum((meep_cs - miepy_cs) ** 2) / meep_norm)
        else:
            rel_L2 = np.nan  # signal too small for meaningful relative error

        results[cs_type] = {
            "rel_L2": rel_L2,
            "max_rel": np.max(rel_err),
            "mean_rel": np.mean(rel_err),
            "meep_cs": meep_cs,
            "miepy_cs": miepy_cs,
            "rel_err": rel_err,
        }

    return results


def print_errors(test_name, errors):
    """Print error summary for a test case."""
    print(f"\n  {'Type':<6} {'Rel L2':>10} {'Max Rel':>10} {'Mean Rel':>10}")
    print(f"  {'-' * 6} {'-' * 10} {'-' * 10} {'-' * 10}")
    for cs_type in ["scat", "abs", "ext"]:
        e = errors[cs_type]
        l2 = f"{e['rel_L2']:10.4f}" if not np.isnan(e["rel_L2"]) else "       N/A"
        print(f"  {cs_type:<6} {l2} {e['max_rel']:10.4f} {e['mean_rel']:10.4f}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_comparison(test_name, meep_data, miepy_data, errors):
    """Plot cross-section comparison (Meep dots vs miepy lines)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    meep_wl = 1.0 / meep_data["frequency"] / nm  # wavelengths in nm
    miepy_wl = miepy_data["wavelengths"] / nm

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_title(test_name, fontsize=16)

    colors = {"scat": "C0", "abs": "C1", "ext": "C2"}
    labels = {"scat": "Scattering", "abs": "Absorption", "ext": "Extinction"}

    # Only plot absorption if signal is meaningful (not just noise)
    meep_ext_max = np.max(np.abs(meep_data["cross_section_ext"]))
    meep_abs_max = np.max(np.abs(meep_data["cross_section_abs"]))
    has_absorption = meep_abs_max > 0.01 * meep_ext_max

    plot_types = ["scat", "abs", "ext"] if has_absorption else ["scat", "ext"]

    for cs_type in plot_types:
        c = colors[cs_type]
        label = labels[cs_type]

        # Meep as circles
        ax.plot(meep_wl, meep_data[f"cross_section_{cs_type}"], "o", color=c, markersize=4, label=f"{label} (Meep)")
        # miepy as lines
        ax.plot(miepy_wl, miepy_data[f"cross_section_{cs_type}"], "-", color=c, label=f"{label} (miepy)")

    ax.set_ylabel("Cross-section (m²)", fontsize=15)
    ax.set_xlabel("Wavelength (nm)", fontsize=15)
    ax.tick_params(labelsize=13)
    ax.legend(fontsize=12, ncol=2)
    ax.axhline(0, ls="--", color="black", alpha=0.3)

    plt.tight_layout()
    out_path = CACHE_DIR / f"{test_name}_comparison.svg"
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  Plot saved to {out_path}")


# ---------------------------------------------------------------------------
# Test case definitions
# ---------------------------------------------------------------------------


def _rotation_matrix(theta, phi):
    """Compute rotation matrix from spherical angles via quaternion."""
    import quaternion as quat

    q = quat.from_spherical_coords(theta, phi)
    return quat.as_rotation_matrix(q), q


def _meep_material_dielectric(index=3.5):
    import meep

    return meep.Medium(index=index)


# --- Dielectric test cases ---


def _run_sphere_dielectric(resolution):
    import meep

    r = 100 * nm
    mat = _meep_material_dielectric()
    geometry = [meep.Sphere(radius=r, center=meep.Vector3(0, 0, 0), material=mat)]
    box = [2 * r] * 3
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_sphere_dielectric():
    import miepy

    return (
        miepy.sphere([0, 0, 0], radius=100 * nm, material=miepy.constant_material(3.5**2)),
        3,
        None,
    )


def _run_spheroid_oblate(resolution):
    import meep

    a, c = 100 * nm, 50 * nm  # xy=100nm, z=50nm
    mat = _meep_material_dielectric()
    geometry = [
        meep.Ellipsoid(
            center=meep.Vector3(0, 0, 0),
            size=meep.Vector3(2 * a, 2 * a, 2 * c),
            material=mat,
        )
    ]
    box = [2 * a] * 3
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_spheroid_oblate():
    import miepy

    return (
        miepy.spheroid([0, 0, 0], axis_xy=100 * nm, axis_z=50 * nm, material=miepy.constant_material(3.5**2)),
        3,
        None,
    )


def _run_spheroid_prolate(resolution):
    import meep

    a, c = 40 * nm, 200 * nm
    mat = _meep_material_dielectric()
    geometry = [
        meep.Ellipsoid(
            center=meep.Vector3(0, 0, 0),
            size=meep.Vector3(2 * a, 2 * a, 2 * c),
            material=mat,
        )
    ]
    box = [2 * c] * 3  # use max dimension
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_spheroid_prolate():
    import miepy

    return (
        miepy.spheroid([0, 0, 0], axis_xy=40 * nm, axis_z=200 * nm, material=miepy.constant_material(3.5**2)),
        4,
        None,
    )


def _run_spheroid_extreme(resolution):
    import meep

    a, c = 30 * nm, 300 * nm
    mat = _meep_material_dielectric()
    geometry = [
        meep.Ellipsoid(
            center=meep.Vector3(0, 0, 0),
            size=meep.Vector3(2 * a, 2 * a, 2 * c),
            material=mat,
        )
    ]
    box = [2 * c] * 3
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_spheroid_extreme():
    import miepy

    return (
        miepy.spheroid(
            [0, 0, 0],
            axis_xy=30 * nm,
            axis_z=300 * nm,
            material=miepy.constant_material(3.5**2),
            extended_precision=True,
        ),
        7,
        None,
    )


def _run_cylinder_aligned(resolution):
    import meep

    r, h = 75 * nm, 150 * nm
    mat = _meep_material_dielectric()
    geometry = [
        meep.Cylinder(
            radius=r,
            height=h,
            axis=meep.Vector3(0, 0, 1),
            center=meep.Vector3(0, 0, 0),
            material=mat,
        )
    ]
    enc_r = np.sqrt(r**2 + (h / 2) ** 2)
    box = [2 * enc_r] * 3
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_cylinder_aligned():
    import miepy

    return (
        miepy.cylinder([0, 0, 0], radius=75 * nm, height=150 * nm, material=miepy.constant_material(3.5**2)),
        3,
        None,
    )


def _run_cylinder_tilted(resolution):
    import meep

    r, h = 75 * nm, 150 * nm
    theta, phi = np.pi / 6, np.pi / 4
    R, q = _rotation_matrix(theta, phi)
    mat = _meep_material_dielectric()

    axis = meep.Vector3(*R[:, 2])
    geometry = [
        meep.Cylinder(
            radius=r,
            height=h,
            axis=axis,
            center=meep.Vector3(0, 0, 0),
            material=mat,
        )
    ]
    enc_r = np.sqrt(r**2 + (h / 2) ** 2)
    box = [2 * enc_r] * 3
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_cylinder_tilted():
    import miepy
    import quaternion as quat

    theta, phi = np.pi / 6, np.pi / 4
    q = quat.from_spherical_coords(theta, phi)
    return (
        miepy.cylinder(
            [0, 0, 0], radius=75 * nm, height=150 * nm, material=miepy.constant_material(3.5**2), orientation=q
        ),
        3,
        None,
    )


def _run_ellipsoid_aligned(resolution):
    import meep

    rx, ry, rz = 100 * nm, 75 * nm, 50 * nm
    mat = _meep_material_dielectric()
    geometry = [
        meep.Ellipsoid(
            center=meep.Vector3(0, 0, 0),
            size=meep.Vector3(2 * rx, 2 * ry, 2 * rz),
            material=mat,
        )
    ]
    enc_r = max(rx, ry, rz)
    box = [2 * enc_r] * 3
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_ellipsoid_aligned():
    import miepy

    return (
        miepy.ellipsoid([0, 0, 0], rx=100 * nm, ry=75 * nm, rz=50 * nm, material=miepy.constant_material(3.5**2)),
        3,
        None,
    )


def _run_ellipsoid_rotated(resolution):
    import meep

    rx, ry, rz = 100 * nm, 75 * nm, 50 * nm
    theta, phi = np.pi / 4, np.pi / 3
    R, q = _rotation_matrix(theta, phi)
    mat = _meep_material_dielectric()

    geometry = [
        meep.Ellipsoid(
            center=meep.Vector3(0, 0, 0),
            size=meep.Vector3(2 * rx, 2 * ry, 2 * rz),
            e1=meep.Vector3(*R[:, 0]),
            e2=meep.Vector3(*R[:, 1]),
            e3=meep.Vector3(*R[:, 2]),
            material=mat,
        )
    ]
    enc_r = max(rx, ry, rz)
    box = [2 * enc_r] * 3
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_ellipsoid_rotated():
    import miepy
    import quaternion as quat

    theta, phi = np.pi / 4, np.pi / 3
    q = quat.from_spherical_coords(theta, phi)
    return (
        miepy.ellipsoid(
            [0, 0, 0], rx=100 * nm, ry=75 * nm, rz=50 * nm, material=miepy.constant_material(3.5**2), orientation=q
        ),
        3,
        None,
    )


def _run_cube_dielectric(resolution):
    import meep

    W = 200 * nm
    mat = _meep_material_dielectric()
    geometry = [
        meep.Block(
            center=meep.Vector3(0, 0, 0),
            size=meep.Vector3(W, W, W),
            material=mat,
        )
    ]
    box = [W] * 3
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm)
    return sim.run()


def _miepy_cube_dielectric():
    import miepy

    return (
        miepy.cube([0, 0, 0], width=200 * nm, material=miepy.constant_material(3.5**2)),
        6,
        None,
    )


def _run_hexprism_dielectric(resolution):
    import meep

    W = 120 * nm  # side length
    H = 150 * nm
    N = 6

    # Compute vertices in um (Meep's Prism has numerical issues at SI-meter scale)
    radius = W / np.sqrt(2 * (1 - np.cos(2 * np.pi / N)))
    radius_um = radius / um
    H_um = H / um
    angles = np.arange(N) * 2 * np.pi / N
    vertices = [meep.Vector3(radius_um * np.cos(a), radius_um * np.sin(a)) for a in angles]

    mat = _meep_material_dielectric()
    geometry = [
        meep.Prism(
            vertices=vertices,
            height=H_um,
            center=meep.Vector3(0, 0, 0),
            material=mat,
        )
    ]
    enc_r = np.sqrt(radius**2 + (H / 2) ** 2)
    box = [2 * enc_r] * 3  # in meters — MeepScatteringSim converts internally
    sim = MeepScatteringSim(geometry, box, resolution, 400 * nm, 1000 * nm, length_scale=um)
    return sim.run()


def _miepy_hexprism_dielectric():
    import miepy

    return (
        miepy.regular_prism([0, 0, 0], N=6, width=120 * nm, height=150 * nm, material=miepy.constant_material(3.5**2)),
        6,
        None,
    )


# --- Gold test cases ---


def _get_gold_eps_data(wavelengths_si):
    """Query meep.materials.Au for permittivity at given SI wavelengths.

    Returns complex eps array.
    """
    import meep
    import meep.materials

    eps_values = np.zeros(len(wavelengths_si), dtype=complex)
    for i, wl in enumerate(wavelengths_si):
        freq_um = um / wl  # frequency in cycles/um (Meep's a=1um convention)
        eps_values[i] = meep.materials.Au.epsilon(freq_um)[0][0]
    return eps_values


def _run_sphere_gold(resolution):
    import meep
    import meep.materials

    r_si = 75 * nm
    r_um = r_si / um

    geometry = [
        meep.Sphere(
            radius=r_um,
            center=meep.Vector3(0, 0, 0),
            material=meep.materials.Au,
        )
    ]
    box_um = [2 * r_um] * 3
    res_um = resolution * um  # convert from pixels/m to pixels/um

    sim = MeepScatteringSimGold(geometry, box_um, res_um, 0.4, 1.0, nfreq=40)
    return sim.run()


def _miepy_sphere_gold():
    """Return (particle, lmax, medium) using eps data from Meep's Au model."""
    import miepy

    wavelengths = np.linspace(400 * nm, 1000 * nm, 100)
    eps = _get_gold_eps_data(wavelengths)
    Au = miepy.data_material(wavelengths, eps)

    return (
        miepy.sphere([0, 0, 0], radius=75 * nm, material=Au),
        3,
        None,
    )


def _run_spheroid_gold(resolution):
    import meep
    import meep.materials

    a_si, c_si = 75 * nm, 50 * nm
    a_um, c_um = a_si / um, c_si / um

    geometry = [
        meep.Ellipsoid(
            center=meep.Vector3(0, 0, 0),
            size=meep.Vector3(2 * a_um, 2 * a_um, 2 * c_um),
            material=meep.materials.Au,
        )
    ]
    box_um = [2 * a_um] * 3
    res_um = resolution * um

    sim = MeepScatteringSimGold(geometry, box_um, res_um, 0.4, 1.0, nfreq=40)
    return sim.run()


def _miepy_spheroid_gold():
    import miepy

    wavelengths = np.linspace(400 * nm, 1000 * nm, 100)
    eps = _get_gold_eps_data(wavelengths)
    Au = miepy.data_material(wavelengths, eps)

    return (
        miepy.spheroid([0, 0, 0], axis_xy=75 * nm, axis_z=50 * nm, material=Au),
        3,
        None,
    )


# ---------------------------------------------------------------------------
# Test case registry
# ---------------------------------------------------------------------------

TEST_CASES = {
    "sphere_dielectric": {
        "description": "Sphere r=100nm, n=3.5 — Mie baseline",
        "meep_func": _run_sphere_dielectric,
        "miepy_func": _miepy_sphere_dielectric,
    },
    "spheroid_oblate": {
        "description": "Oblate spheroid 2:1, xy=100nm z=50nm, n=3.5",
        "meep_func": _run_spheroid_oblate,
        "miepy_func": _miepy_spheroid_oblate,
    },
    "spheroid_prolate": {
        "description": "Prolate spheroid 5:1, xy=40nm z=200nm, n=3.5",
        "meep_func": _run_spheroid_prolate,
        "miepy_func": _miepy_spheroid_prolate,
    },
    "spheroid_extreme": {
        "description": "Prolate spheroid 10:1, xy=30nm z=300nm, n=3.5",
        "meep_func": _run_spheroid_extreme,
        "miepy_func": _miepy_spheroid_extreme,
    },
    "cylinder_aligned": {
        "description": "Cylinder r=75nm h=150nm, n=3.5, aligned",
        "meep_func": _run_cylinder_aligned,
        "miepy_func": _miepy_cylinder_aligned,
    },
    "cylinder_tilted": {
        "description": "Cylinder r=75nm h=150nm, n=3.5, tilted theta=pi/6 phi=pi/4",
        "meep_func": _run_cylinder_tilted,
        "miepy_func": _miepy_cylinder_tilted,
    },
    "ellipsoid_aligned": {
        "description": "Ellipsoid 100/75/50nm, n=3.5, aligned",
        "meep_func": _run_ellipsoid_aligned,
        "miepy_func": _miepy_ellipsoid_aligned,
    },
    "ellipsoid_rotated": {
        "description": "Ellipsoid 100/75/50nm, n=3.5, rotated theta=pi/4 phi=pi/3",
        "meep_func": _run_ellipsoid_rotated,
        "miepy_func": _miepy_ellipsoid_rotated,
    },
    "cube_dielectric": {
        "description": "Cube W=200nm, n=3.5",
        "meep_func": _run_cube_dielectric,
        "miepy_func": _miepy_cube_dielectric,
    },
    "hexprism_dielectric": {
        "description": "Hexagonal prism W=120nm H=150nm, n=3.5",
        "meep_func": _run_hexprism_dielectric,
        "miepy_func": _miepy_hexprism_dielectric,
    },
    "sphere_gold": {
        "description": "Gold sphere r=75nm, Drude-Lorentz",
        "meep_func": _run_sphere_gold,
        "miepy_func": _miepy_sphere_gold,
    },
    "spheroid_gold": {
        "description": "Gold oblate spheroid xy=75nm z=50nm, Drude-Lorentz",
        "meep_func": _run_spheroid_gold,
        "miepy_func": _miepy_spheroid_gold,
    },
}


# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------


def run_test(test_name, resolution, force=False, do_plot=False):
    """Run a single test case: Meep + miepy comparison."""
    tc = TEST_CASES[test_name]
    print(f"\n{'=' * 60}")
    print(f"Test: {test_name}")
    print(f"  {tc['description']}")
    res_nm = round(1 / (resolution * nm))
    print(f"  Resolution: {res_nm} nm grid spacing")

    # Run or load Meep
    meep_data = load_or_run(test_name, resolution, tc["meep_func"], force=force)

    # Run miepy
    print(f"  Running miepy...")
    particle, lmax, medium = tc["miepy_func"]()
    wavelengths = np.linspace(400 * nm, 1000 * nm, 100)
    C_scat, C_abs, C_ext = compute_miepy_cross_sections(particle, wavelengths, lmax, medium=medium)
    miepy_data = {
        "wavelengths": wavelengths,
        "cross_section_scat": C_scat,
        "cross_section_abs": C_abs,
        "cross_section_ext": C_ext,
    }

    # Compare
    errors = compute_errors(meep_data, miepy_data, meep_data["frequency"])
    print_errors(test_name, errors)

    if do_plot:
        plot_comparison(test_name, meep_data, miepy_data, errors)

    return errors


def run_summary(do_plot=False):
    """Print summary table of all cached results."""
    print(f"\n{'=' * 70}")
    print(f"FDTD Validation Summary")
    print(f"{'=' * 70}")

    # Find all cached files
    if not CACHE_DIR.exists():
        print("No cached results found.")
        return

    npz_files = sorted(CACHE_DIR.glob("*.npz"))
    if not npz_files:
        print("No cached results found.")
        return

    # Group by test name, preferring finest resolution (smallest nm value)
    cached_tests = {}
    cached_res = {}
    for f in npz_files:
        name = f.stem
        # Parse test_name_res{N}nm
        parts = name.rsplit("_res", 1)
        if len(parts) == 2:
            test_name = parts[0]
            try:
                res_nm = int(parts[1].replace("nm", ""))
            except ValueError:
                continue
            if test_name not in cached_tests or res_nm < cached_res[test_name]:
                cached_tests[test_name] = f
                cached_res[test_name] = res_nm

    if not cached_tests:
        print("No valid cached results found.")
        return

    print(f"\n{'Test':<25} {'Res':>5} {'Scat L2':>10} {'Abs L2':>10} {'Ext L2':>10}")
    print(f"{'-' * 25} {'-' * 5} {'-' * 10} {'-' * 10} {'-' * 10}")

    for test_name in TEST_CASES:
        if test_name not in cached_tests:
            continue

        f = cached_tests[test_name]
        meep_data = {k: v for k, v in np.load(f).items()}

        # Parse resolution from filename
        res_str = f.stem.rsplit("_res", 1)[1].replace("nm", "")
        res_nm = int(res_str)

        # Run miepy for comparison
        tc = TEST_CASES[test_name]
        particle, lmax, medium = tc["miepy_func"]()
        wavelengths = np.linspace(400 * nm, 1000 * nm, 100)
        C_scat, C_abs, C_ext = compute_miepy_cross_sections(particle, wavelengths, lmax, medium=medium)
        miepy_data = {
            "wavelengths": wavelengths,
            "cross_section_scat": C_scat,
            "cross_section_abs": C_abs,
            "cross_section_ext": C_ext,
        }

        errors = compute_errors(meep_data, miepy_data, meep_data["frequency"])

        if do_plot:
            plot_comparison(test_name, meep_data, miepy_data, errors)

        def _fmt(v):
            return f"{v:10.4f}" if not np.isnan(v) else "       N/A"

        print(
            f"{test_name:<25} {res_nm:>4}nm "
            f"{_fmt(errors['scat']['rel_L2'])} "
            f"{_fmt(errors['abs']['rel_L2'])} "
            f"{_fmt(errors['ext']['rel_L2'])}"
        )


def main():
    parser = argparse.ArgumentParser(description="FDTD (Meep) validation for miepy T-matrix solver")
    parser.add_argument("--list", action="store_true", help="List all test cases")
    parser.add_argument("--run", nargs="+", metavar="TEST", help="Run specified test cases (or 'all')")
    parser.add_argument(
        "--resolution",
        default="standard",
        choices=RESOLUTION_PRESETS.keys(),
        help="Resolution preset (default: standard)",
    )
    parser.add_argument("--plot", action="store_true", help="Generate comparison plots")
    parser.add_argument("--force", action="store_true", help="Ignore cached results and re-run Meep")
    parser.add_argument("--summary", action="store_true", help="Show summary of all cached results")

    args = parser.parse_args()

    if args.list:
        print("Available test cases:")
        print(f"  {'Name':<25} Description")
        print(f"  {'-' * 25} {'-' * 50}")
        for name, tc in TEST_CASES.items():
            print(f"  {name:<25} {tc['description']}")
        return

    if args.summary:
        run_summary(do_plot=args.plot)
        return

    if args.run:
        resolution = RESOLUTION_PRESETS[args.resolution]

        test_names = args.run
        if "all" in test_names:
            test_names = list(TEST_CASES.keys())

        for name in test_names:
            if name not in TEST_CASES:
                print(f"Unknown test case: {name}")
                print(f"Use --list to see available test cases.")
                sys.exit(1)

        for name in test_names:
            run_test(name, resolution, force=args.force, do_plot=args.plot)

        return

    parser.print_help()


if __name__ == "__main__":
    main()
