from get_tmatrix import nfmds_solver

def tmatrix_sphere():
    pass

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
                  index=eps**0.5, index_m=eps_m**0.5, complex_plane=complex_plane)
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
                  index=eps**0.5, index_m=eps_m**0.5, complex_plane=complex_plane)
    parameters.update(kwargs)

    return nfmds_solver(lmax, parameters, extended_precision=extended_precision)
