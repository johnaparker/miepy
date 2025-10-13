"""MiePy.
=======
Python module to calcuate scattering coefficients of a plane wave incident on a sphere or core-shell structure using Mie theory
"""

# main submodules
import warnings

from . import (
    constants,
    coordinates,
    cpp,
    flux,
    forces,
    interactions,
    material_functions,
    microscope,
    particles,
    sources,
    symmetry,
    tmatrix,
    utils,
    vsh,
)
from .cluster import cluster
from .cpp.interactions import solver
from .interface import interface
from .material_functions.create import constant_material, data_material, dielectric, function_material
from .materials.predefined import materials
from .microscope import cluster_microscope, microscope
from .mie_single.mie_core_shell import single_mie_core_shell
from .mie_single.mie_sphere import single_mie_sphere
from .mie_single.scattering import (
    absorbption_per_multipole,
    cross_sections,
    extinction_per_multipole,
    multipole_label,
    scattering_per_multipole,
)
from .particles import core_shell, cube, cylinder, ellipsoid, regular_prism, sphere, sphere_cluster_particle, spheroid
from .sphere_cluster import sphere_cluster
from .visual.view3d import visualize
from .vsh import VSH, cluster_coefficients, expand_E, expand_E_far, expand_H, expand_H_far, mode_indices, vsh_mode
from .vsh.mode_indices import reduced_index

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import quaternion
