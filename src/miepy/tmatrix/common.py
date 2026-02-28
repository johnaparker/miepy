import numpy as np

import miepy

from .get_tmatrix import nfmds_solver, tmatrix_solvers


def tmatrix_sphere(radius, wavelength, eps, eps_m, lmax, conducting=False):
    """Compute the T-matrix of a sphere, using regular Mie theory.

    Arguments:
        radius      sphere radius
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        conducting  if True, calculate for conducting sphere (default: False)
    """
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    tmatrix = np.zeros([2, rmax, 2, rmax], dtype=complex)
    k_medium = 2 * np.pi * eps_m**0.5 / wavelength

    for i, n, _m in miepy.mode_indices(lmax):
        an, bn = miepy.mie_single.mie_sphere_scattering_coefficients(
            radius, n, eps, 1, eps_m, 1, k_medium, conducting=conducting
        )
        tmatrix[0, i, 0, i] = an
        tmatrix[1, i, 1, i] = bn

    return tmatrix


def tmatrix_core_shell(radius, thickness, wavelength, eps_core, eps_shell, eps_m, lmax):
    """Compute the T-matrix of a core-shell, using regular Mie theory.

    Arguments:
        radius      core radius
        wavelength  incident wavelength
        eps_core    particle permittivity
        eps_shell  shell permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
    """
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    tmatrix = np.zeros([2, rmax, 2, rmax], dtype=complex)
    particle = miepy.single_mie_core_shell(
        radius,
        radius + thickness,
        material_in=miepy.dielectric(eps=eps_core),
        material_out=miepy.dielectric(eps=eps_shell),
        medium=miepy.dielectric(eps=eps_m),
        lmax=lmax,
        wavelength=wavelength,
    )

    particle.solve()

    for i, n, _m in miepy.mode_indices(lmax):
        tmatrix[0, i, 0, i] = particle.an[0, n - 1]
        tmatrix[1, i, 1, i] = particle.bn[0, n - 1]

    return tmatrix


def tmatrix_spheroid(axis_xy, axis_z, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of a spheroid.

    Arguments:
        axis_xy     length of semiaxes perpendicular to the axis of symmetry
        axis_z      length of semiaxis along axis of symmetry
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to C++ or Fortran solver
    """
    conducting = kwargs.get('conducting', False)
    k = 2 * np.pi * eps_m**0.5 / wavelength
    if conducting:
        n_rel = complex(1.0 / eps_m**0.5)
    else:
        n_rel = complex((eps / eps_m) ** 0.5)
    complex_plane = axis_xy > axis_z
    Nint = kwargs.get('Nint', 200)
    eps_z = kwargs.get('eps_z_re_im', 0.95)
    use_ds = kwargs.get('use_ds', True)

    return miepy.cpp.tmatrix.compute_spheroid(
        axis_z, axis_xy, k, n_rel, lmax, Nint,
        use_ds, complex_plane, eps_z, extended_precision, conducting)


def tmatrix_cylinder(radius, height, wavelength, eps, eps_m, lmax, rounded=False, extended_precision=False, **kwargs):
    """Compute the T-matrix of a cylinder, with sharp or rounded (if oblate) edges.

    Arguments:
        radius      radius of cylinder
        height      height of cylinder
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        rounded (bool)    if True, and cylinder is oblate, the cylinder's edges are rounded (default: False)
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to C++ or Fortran solver
    """
    if height >= 2 * radius and rounded:
        raise ValueError("prolate cylinders (height >= diameter) cannot be rounded")

    conducting = kwargs.get('conducting', False)
    k = 2 * np.pi * eps_m**0.5 / wavelength
    if conducting:
        n_rel = complex(1.0 / eps_m**0.5)
    else:
        n_rel = complex((eps / eps_m) ** 0.5)
    complex_plane = 2 * radius > height
    Nint = kwargs.get('Nint', 200)
    eps_z = kwargs.get('eps_z_re_im', 0.95)
    use_ds = kwargs.get('use_ds', True)

    if rounded:
        return miepy.cpp.tmatrix.compute_rounded_cylinder(
            height / 2, radius, k, n_rel, lmax, Nint,
            use_ds, complex_plane, eps_z, extended_precision, conducting)
    else:
        return miepy.cpp.tmatrix.compute_cylinder(
            height / 2, radius, k, n_rel, lmax, Nint,
            use_ds, complex_plane, eps_z, extended_precision, conducting)


def tmatrix_ellipsoid(rx, ry, rz, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of an ellipsoid.

    Arguments:
        rx,ry,rz    radii of the 3 axes
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to C++ solver
    """
    conducting = kwargs.get('conducting', False)
    k = 2 * np.pi * eps_m**0.5 / wavelength
    if conducting:
        n_rel = complex(1.0 / eps_m**0.5)
    else:
        n_rel = complex((eps / eps_m) ** 0.5)
    Nint1 = kwargs.get('Nint1', 50)
    Nint2 = kwargs.get('Nint2', 50)

    return miepy.cpp.tmatrix.compute_ellipsoid(
        rx, ry, rz, k, n_rel, lmax, Nint1, Nint2,
        extended_precision, conducting)


def tmatrix_square_prism(side, height, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of a square prism (cube when side==height).

    Arguments:
        side        side width of the prism
        height      height of the prism
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to C++ solver
    """
    return tmatrix_regular_prism(4, side, height, wavelength, eps, eps_m, lmax,
                                 extended_precision=extended_precision, **kwargs)


def tmatrix_regular_prism(N, side, height, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of a regular N-hedral prism.

    Arguments:
        N           number of vertices
        side        side width of the prism
        height      height of the prism
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to C++ solver
    """
    conducting = kwargs.get('conducting', False)
    k = 2 * np.pi * eps_m**0.5 / wavelength
    if conducting:
        n_rel = complex(1.0 / eps_m**0.5)
    else:
        n_rel = complex((eps / eps_m) ** 0.5)
    Nint1 = kwargs.get('Nint1', 50)
    Nint2 = kwargs.get('Nint2', 50)

    return miepy.cpp.tmatrix.compute_regular_prism(
        N, side / 2, height / 2, k, n_rel, lmax, Nint1, Nint2,
        extended_precision, conducting)


def tmatrix_sphere_cluster(pos, radii, lmax, lmax_cluster, wavelength, eps, eps_m, extended_precision=False, **kwargs):
    parameters = dict(
        pos=pos, radii=radii, Nrank_particles=lmax, wavelength=wavelength, index=eps**0.5, index_m=eps_m**0.5
    )
    parameters.update(kwargs)

    return nfmds_solver(
        lmax_cluster, parameters, solver=tmatrix_solvers.sphere_cluster, extended_precision=extended_precision
    )
