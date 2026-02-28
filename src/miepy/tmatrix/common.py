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


def _build_source_matrix(positions, origin, lmax_p, lmax_out, k):
    """Build source matrix translating origin VSH modes to particle positions.

    For each origin VSH mode (n,m) with polarization a, computes the incident
    field at each particle position using the VSH addition theorem.

    Arguments:
        positions[N,3]  particle positions
        origin[3]       cluster origin
        lmax_p          internal lmax for particles
        lmax_out        lmax for origin modes
        k               medium wavenumber

    Returns:
        S[N*2*rmax_p, 2*rmax_out]
    """
    from miepy.cpp.vsh_translation import vsh_translation as _vsh_translation

    Nparticles = positions.shape[0]
    rmax_p = miepy.vsh.lmax_to_rmax(lmax_p)
    rmax_out = miepy.vsh.lmax_to_rmax(lmax_out)

    S = np.zeros((Nparticles * 2 * rmax_p, 2 * rmax_out), dtype=complex)

    for i in range(Nparticles):
        d = positions[i] - origin
        dist = np.linalg.norm(d)

        if dist < 1e-15:
            lmax_min = min(lmax_p, lmax_out)
            rmax_min = miepy.vsh.lmax_to_rmax(lmax_min)
            for a in range(2):
                for r in range(rmax_min):
                    S[i * 2 * rmax_p + a * rmax_p + r, a * rmax_out + r] = 1.0
            continue

        rad, theta, phi = miepy.coordinates.cart_to_sph(*d)

        for s, v, u in miepy.mode_indices(lmax_p):
            for r, n, m in miepy.mode_indices(lmax_out):
                A, B = _vsh_translation(u, v, m, n, rad, theta, phi, k, miepy.vsh_mode.incident)

                S[i * 2 * rmax_p + 0 * rmax_p + s, 0 * rmax_out + r] += A
                S[i * 2 * rmax_p + 0 * rmax_p + s, 1 * rmax_out + r] += B
                S[i * 2 * rmax_p + 1 * rmax_p + s, 0 * rmax_out + r] += B
                S[i * 2 * rmax_p + 1 * rmax_p + s, 1 * rmax_out + r] += A

    return S


def _build_reexpansion_matrix(positions, origin, lmax_p, lmax_out, k):
    """Build re-expansion matrix translating particle scattered modes to origin.

    Same structure as cluster_coefficients() but stored as an explicit matrix.

    Arguments:
        positions[N,3]  particle positions
        origin[3]       cluster origin
        lmax_p          internal lmax for particles
        lmax_out        lmax for origin modes
        k               medium wavenumber

    Returns:
        R[2*rmax_out, N*2*rmax_p]
    """
    from miepy.cpp.vsh_translation import vsh_translation as _vsh_translation

    Nparticles = positions.shape[0]
    rmax_p = miepy.vsh.lmax_to_rmax(lmax_p)
    rmax_out = miepy.vsh.lmax_to_rmax(lmax_out)

    R = np.zeros((2 * rmax_out, Nparticles * 2 * rmax_p), dtype=complex)

    for i in range(Nparticles):
        d = origin - positions[i]
        dist = np.linalg.norm(d)

        if dist < 1e-15:
            lmax_min = min(lmax_p, lmax_out)
            rmax_min = miepy.vsh.lmax_to_rmax(lmax_min)
            for a in range(2):
                for r in range(rmax_min):
                    R[a * rmax_out + r, i * 2 * rmax_p + a * rmax_p + r] = 1.0
            continue

        rad, theta, phi = miepy.coordinates.cart_to_sph(*d)

        for r, n, m in miepy.mode_indices(lmax_out):
            for s, v, u in miepy.mode_indices(lmax_p):
                A, B = _vsh_translation(m, n, u, v, rad, theta, phi, k, miepy.vsh_mode.incident)

                R[0 * rmax_out + r, i * 2 * rmax_p + 0 * rmax_p + s] += A
                R[0 * rmax_out + r, i * 2 * rmax_p + 1 * rmax_p + s] += B
                R[1 * rmax_out + r, i * 2 * rmax_p + 0 * rmax_p + s] += B
                R[1 * rmax_out + r, i * 2 * rmax_p + 1 * rmax_p + s] += A

    return R


def tmatrix_sphere_cluster(pos, radii, lmax, lmax_cluster, wavelength, eps, eps_m, **kwargs):
    """Compute the T-matrix of a rigid sphere cluster.

    Uses Mie theory, VSH translations, and the aggregate T-matrix formulation
    to compute the cluster T-matrix without external Fortran code.

    The T-matrix maps incident VSH modes at the cluster origin to scattered
    VSH modes: T = R @ M @ solve(I + T_agg, S)

    Arguments:
        pos[N,3]        positions of the spheres
        radii[N]        radii of the spheres (or scalar)
        lmax            per-particle lmax values (scalar or array)
        lmax_cluster    lmax for the output T-matrix
        wavelength      incident wavelength
        eps[N]          particle permittivities (or scalar)
        eps_m           medium permittivity
        kwargs          additional keywords (ignored)
    """
    pos = np.asarray(pos, dtype=float)
    radii = np.atleast_1d(np.asarray(radii, dtype=float))
    lmax_arr = np.atleast_1d(np.asarray(lmax, dtype=int))
    eps_arr = np.atleast_1d(np.asarray(eps, dtype=complex))

    Nparticles = pos.shape[0]
    if len(radii) == 1:
        radii = np.full(Nparticles, radii[0])
    if len(lmax_arr) == 1:
        lmax_arr = np.full(Nparticles, lmax_arr[0], dtype=int)
    if len(eps_arr) == 1:
        eps_arr = np.full(Nparticles, eps_arr[0])

    lmax_internal = int(np.max(lmax_arr))
    k = 2 * np.pi * eps_m**0.5 / wavelength
    origin = np.zeros(3)
    rmax_p = miepy.vsh.lmax_to_rmax(lmax_internal)
    rmax_out = miepy.vsh.lmax_to_rmax(lmax_cluster)

    # Step 1: Compute Mie scattering coefficients
    mie = np.zeros([Nparticles, 2, lmax_internal], dtype=complex)
    for i in range(Nparticles):
        for n in range(1, int(lmax_arr[i]) + 1):
            an, bn = miepy.mie_single.mie_sphere_scattering_coefficients(
                radii[i], n, eps_arr[i], 1, eps_m, 1, k
            )
            mie[i, 0, n - 1] = an
            mie[i, 1, n - 1] = bn

    # Step 2: Build aggregate T-matrix and form (I + T_agg)
    T_agg = miepy.interactions.sphere_aggregate_tmatrix(pos, mie, k)
    size = Nparticles * 2 * rmax_p
    A_matrix = np.eye(size, dtype=complex) + T_agg.reshape(size, size)

    # Step 3: Build source matrix (origin → particles)
    S = _build_source_matrix(pos, origin, lmax_internal, lmax_cluster, k)

    # Step 4: Solve (I + T_agg) @ X = S for all RHS at once
    X = np.linalg.solve(A_matrix, S)

    # Step 5: Apply Mie scattering (diagonal operation)
    mie_diag = np.zeros(size, dtype=complex)
    for i in range(Nparticles):
        for a in range(2):
            for r, n, _m in miepy.mode_indices(lmax_internal):
                mie_diag[i * 2 * rmax_p + a * rmax_p + r] = mie[i, a, n - 1]
    Y = mie_diag[:, np.newaxis] * X

    # Step 6: Re-expand to origin
    R = _build_reexpansion_matrix(pos, origin, lmax_internal, lmax_cluster, k)
    T_flat = R @ Y

    return T_flat.reshape(2, rmax_out, 2, rmax_out)
