"""Tests for core-shell T-matrix consistency with single Mie theory."""

import numpy as np

import miepy

nm = 1e-9

# wavelength from 400nm to 1000nm
wavelengths = np.linspace(400 * nm, 1000 * nm, 10)

core_radius = 50 * nm
shell_thickness = 30 * nm
lmax = 3

# Ag shell on dielectric core — absorbing shell exercises absorption cross-section
core_material = miepy.constant_material(2.5**2)
shell_material = miepy.materials.Ag()
medium = miepy.materials.water()

# Single Mie core-shell theory (reference)
mie = miepy.single_mie_core_shell(
    core_radius,
    core_radius + shell_thickness,
    material_in=core_material,
    material_out=shell_material,
    wavelength=wavelengths,
    lmax=lmax,
    medium=medium,
)
S_ref, A_ref, E_ref = mie.cross_sections()

# GMT (cluster) with a single core-shell particle
source = miepy.sources.plane_wave.from_string(polarization="x")

scat = np.zeros_like(wavelengths)
absorb = np.zeros_like(wavelengths)
extinct = np.zeros_like(wavelengths)

for i, wavelength in enumerate(wavelengths):
    system = miepy.cluster(
        particles=[
            miepy.core_shell(
                position=[0, 0, 0],
                core_radius=core_radius,
                shell_thickness=shell_thickness,
                core_material=core_material,
                shell_material=shell_material,
            )
        ],
        source=source,
        wavelength=wavelength,
        lmax=lmax,
        medium=medium,
    )

    scat[i], absorb[i], extinct[i] = system.cross_sections()


def test_core_shell_cross_sections():
    """Cross-sections from cluster with one core-shell must match single Mie theory."""
    for name, cluster_val, ref_val in [
        ("scattering", scat, S_ref),
        ("absorption", absorb, A_ref),
        ("extinction", extinct, E_ref),
    ]:
        L2 = np.linalg.norm(cluster_val - ref_val) / len(wavelengths)
        avg = np.average(np.abs(ref_val) + np.abs(cluster_val)) / 2
        assert L2 < 1e-13 * avg, f"{name} mismatch: L2={L2}, avg={avg}"


def test_core_shell_absorption_nonnegative():
    """Cluster absorption cross-section must be non-negative."""
    assert np.all(absorb >= 0), f"Negative absorption found: {absorb}"


def test_core_shell_degenerate_to_sphere():
    """When core and shell have the same material, T-matrix must match tmatrix_sphere."""
    wavelength = 600 * nm
    outer_radius = core_radius + shell_thickness
    mat = miepy.constant_material(3.0**2)
    eps = mat.eps(wavelength)
    eps_m = medium.eps(wavelength)

    T_cs = miepy.tmatrix.tmatrix_core_shell(
        core_radius, shell_thickness, wavelength, eps, eps, eps_m, lmax
    )
    T_sp = miepy.tmatrix.tmatrix_sphere(outer_radius, wavelength, eps, eps_m, lmax)

    assert np.allclose(T_cs, T_sp, atol=1e-15, rtol=1e-10), (
        f"Core-shell with uniform material differs from sphere:\n"
        f"max diff = {np.max(np.abs(T_cs - T_sp))}"
    )


def test_degenerate_core_shell_single_particle():
    """A single core-shell with uniform material must give the same cross-sections as a sphere."""
    uniform_mat = miepy.constant_material(3.7**2 + 0.1j)
    outer_radius = core_radius + shell_thickness
    test_wavelengths = np.linspace(400 * nm, 800 * nm, 5)

    cs_scat = np.zeros_like(test_wavelengths)
    cs_abs = np.zeros_like(test_wavelengths)
    cs_ext = np.zeros_like(test_wavelengths)
    sp_scat = np.zeros_like(test_wavelengths)
    sp_abs = np.zeros_like(test_wavelengths)
    sp_ext = np.zeros_like(test_wavelengths)

    for i, wl in enumerate(test_wavelengths):
        cs_system = miepy.cluster(
            particles=[
                miepy.core_shell(
                    position=[0, 0, 0],
                    core_radius=core_radius,
                    shell_thickness=shell_thickness,
                    core_material=uniform_mat,
                    shell_material=uniform_mat,
                )
            ],
            source=source,
            wavelength=wl,
            lmax=lmax,
            medium=medium,
        )
        cs_scat[i], cs_abs[i], cs_ext[i] = cs_system.cross_sections()

        sp_system = miepy.sphere_cluster(
            position=[0, 0, 0],
            radius=outer_radius,
            material=uniform_mat,
            source=source,
            wavelength=wl,
            lmax=lmax,
            medium=medium,
        )
        sp_scat[i], sp_abs[i], sp_ext[i] = sp_system.cross_sections()

    for name, cs_val, sp_val in [
        ("scattering", cs_scat, sp_scat),
        ("absorption", cs_abs, sp_abs),
        ("extinction", cs_ext, sp_ext),
    ]:
        L2 = np.linalg.norm(cs_val - sp_val) / len(test_wavelengths)
        avg = np.average(np.abs(cs_val) + np.abs(sp_val)) / 2
        assert L2 < 1e-13 * avg, f"{name} mismatch: L2={L2}, avg={avg}"


def test_degenerate_core_shell_interacting_cluster():
    """A cluster of interacting core-shells with uniform material must match a sphere_cluster."""
    uniform_mat = miepy.constant_material(3.7**2 + 0.1j)
    outer_radius = core_radius + shell_thickness
    sep = 300 * nm
    positions = [[sep / 2, 0, 0], [-sep / 2, 0, 0], [0, sep / 2, 0]]
    test_wavelengths = np.linspace(400 * nm, 800 * nm, 5)

    cs_scat = np.zeros_like(test_wavelengths)
    cs_abs = np.zeros_like(test_wavelengths)
    cs_ext = np.zeros_like(test_wavelengths)
    sp_scat = np.zeros_like(test_wavelengths)
    sp_abs = np.zeros_like(test_wavelengths)
    sp_ext = np.zeros_like(test_wavelengths)

    for i, wl in enumerate(test_wavelengths):
        cs_system = miepy.cluster(
            particles=[
                miepy.core_shell(
                    position=pos,
                    core_radius=core_radius,
                    shell_thickness=shell_thickness,
                    core_material=uniform_mat,
                    shell_material=uniform_mat,
                )
                for pos in positions
            ],
            source=source,
            wavelength=wl,
            lmax=lmax,
            medium=medium,
        )
        cs_scat[i], cs_abs[i], cs_ext[i] = cs_system.cross_sections()

        sp_system = miepy.sphere_cluster(
            position=positions,
            radius=outer_radius,
            material=uniform_mat,
            source=source,
            wavelength=wl,
            lmax=lmax,
            medium=medium,
        )
        sp_scat[i], sp_abs[i], sp_ext[i] = sp_system.cross_sections()

    for name, cs_val, sp_val in [
        ("scattering", cs_scat, sp_scat),
        ("absorption", cs_abs, sp_abs),
        ("extinction", cs_ext, sp_ext),
    ]:
        L2 = np.linalg.norm(cs_val - sp_val) / len(test_wavelengths)
        avg = np.average(np.abs(cs_val) + np.abs(sp_val)) / 2
        assert L2 < 1e-13 * avg, f"{name} mismatch: L2={L2}, avg={avg}"
