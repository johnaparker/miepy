import miepy
import numpy as np
from .get_tmatrix import nfmds_solver, tmatrix_solvers

def tmatrix_sphere(radius, wavelength, eps, eps_m, lmax, conducting=False):
    """Compute the T-matrix of a sphere, using regular Mie theory

    Arguments:
        radius      sphere radius
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        conducting  if True, calculate for conducting sphere (default: False)
    """
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    tmatrix = np.zeros([2,rmax,2,rmax], dtype=complex)
    k_medium = 2*np.pi*eps_m**0.5/wavelength

    for i, n, m in miepy.mode_indices(lmax):
        an, bn = miepy.mie_single.mie_sphere_scattering_coefficients(radius,
                          n, eps, 1, eps_m, 1, k_medium, conducting=conducting)
        tmatrix[0,i,0,i] = an
        tmatrix[1,i,1,i] = bn

    return tmatrix

def tmatrix_spheroid(axis_xy, axis_z, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of a spheroid
    
    Arguments:
        axis_xy     length of semiaxes perpendicular to the axis of symmetry
        axis_z      length of semiaxis along axis of symmetry
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to axisymmetric_file function
    """
    complex_plane = True if axis_xy > axis_z else False
    parameters = dict(geometry_type=1, geometry_parameters=[axis_z, axis_xy], wavelength=wavelength,
                  index=eps**.5, index_m=eps_m**0.5, complex_plane=complex_plane, Nparam=1)
    parameters.update(kwargs)

    return nfmds_solver(lmax, parameters, extended_precision=extended_precision)

def tmatrix_cylinder(radius, height, wavelength, eps, eps_m, lmax, rounded=False, extended_precision=False, **kwargs):
    """Compute the T-matrix of a cylinder, with sharp or rounded (if oblate) edges
    
    Arguments:
        radius      radius of cylinder
        height      height of cylinder
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        rounded (bool)    if True, and cylinder is oblate, the cylinder's edges are rounded (default: False)
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to axisymmetric_file function
    """
    complex_plane = True if 2*radius > height else False
    geometry_type = 3 if rounded else 2
    if height >= 2*radius and rounded:
        raise ValueError('prolate cylinders (height >= diameter) cannot be rounded')

    parameters = dict(geometry_type=geometry_type, geometry_parameters=[height/2, radius], wavelength=wavelength,
                  index=eps**0.5, index_m=eps_m**0.5, complex_plane=complex_plane, Nparam=3)
    parameters.update(kwargs)

    return nfmds_solver(lmax, parameters, extended_precision=extended_precision)

def tmatrix_ellipsoid(rx, ry, rz, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of a spheroid
    
    Arguments:
        rx,ry,rz    radii of the 3 axes
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to axisymmetric_file function
    """
    parameters = dict(geometry_type=1, geometry_parameters=[rx, ry, rz], wavelength=wavelength,
                  index=eps**0.5, index_m=eps_m**0.5, Nparam=1, Mrank=lmax, R_symmetry=0)
    parameters.update(kwargs)

    return nfmds_solver(lmax, parameters, solver=tmatrix_solvers.non_axisymmetric,
                        extended_precision=extended_precision)

def tmatrix_ellipsoid(rx, ry, rz, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of a spheroid
    
    Arguments:
        rx,ry,rz    radii of the 3 axes
        wavelength  incident wavelength
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to axisymmetric_file function
    """
    parameters = dict(geometry_type=1, geometry_parameters=[rx, ry, rz], wavelength=wavelength,
                  index=eps**0.5, index_m=eps_m**0.5, Nparam=1, Mrank=lmax, R_symmetry=0)
    parameters.update(kwargs)

    return nfmds_solver(lmax, parameters, solver=tmatrix_solvers.non_axisymmetric,
                        extended_precision=extended_precision)

def tmatrix_square_prism(side, height, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of a spheroid
    
    Arguments:
        width       side width of the prism
        height      height of the prism
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to axisymmetric_file function
    """
    parameters = dict(geometry_type=2, geometry_parameters=[side/2, height/2], wavelength=wavelength,
                  index=eps**0.5, index_m=eps_m**0.5, Nparam=6, Mrank=lmax, R_symmetry=0)
    parameters.update(kwargs)

    return nfmds_solver(lmax, parameters, solver=tmatrix_solvers.non_axisymmetric,
                        extended_precision=extended_precision)

def tmatrix_regular_prism(N, side, height, wavelength, eps, eps_m, lmax, extended_precision=False, **kwargs):
    """Compute the T-matrix of a spheroid
    
    Arguments:
        N           number of vertices
        width       side width of the prism
        height      height of the prism
        eps         particle permittivity
        eps_m       medium permittivity
        lmax        maximum number of multipoles
        extended_precision (bool)    whether to use extended precision (default: False)
        kwargs      additional keywords passed to axisymmetric_file function
    """
    parameters = dict(geometry_type=3, geometry_parameters=[side/2, height/2], wavelength=wavelength,
                  index=eps**0.5, index_m=eps_m**0.5, Nparam=2, Mrank=lmax, R_symmetry=N)
    parameters.update(kwargs)

    return nfmds_solver(lmax, parameters, solver=tmatrix_solvers.non_axisymmetric,
                        extended_precision=extended_precision)
