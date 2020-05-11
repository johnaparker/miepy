"""
MiePy
=======
Python module to calcuate scattering coefficients of a plane wave incident on a sphere or core-shell structure using Mie theory
"""

# main submodules
from . import cpp
from . import vsh
from . import sources
from . import material_functions
from . import interactions
from . import tmatrix
from . import particles
from . import forces
from . import flux
from . import coordinates
from . import utils
from . import symmetry
from . import constants
from . import microscope

from .material_functions.create import dielectric, constant_material, function_material, data_material
from .materials.predefined import materials
from .cluster import cluster
from .sphere_cluster import sphere_cluster
from .mie_single.scattering import (scattering_per_multipole, absorbption_per_multipole,
                        extinction_per_multipole, cross_sections, multipole_label)
from .mie_single.mie_sphere import single_mie_sphere
from .mie_single.mie_core_shell import single_mie_core_shell
from .vsh import (mode_indices, vsh_mode, VSH, expand_E, expand_E_far, expand_H, expand_H_far,
                  cluster_coefficients)
from .particles import (sphere, spheroid, cylinder, ellipsoid, regular_prism, cube,
                       sphere_cluster_particle)
from .cpp.interactions import solver
from .interface import interface
from .visual.view3d import visualize
from .microscope import microscope, cluster_microscope

from .vsh.mode_indices import reduced_index

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import quaternion
