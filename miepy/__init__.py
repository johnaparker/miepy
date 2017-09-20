"""
MiePy
=======
Python module to calcuate scattering coefficients of a plane wave incident on a sphere or core-shell structure using Mie theory
"""

#main submodules
from . import sources
from . import material_functions
from . import scattering
from . import materials

from .mie_sphere import single_mie_sphere
from .mie_core_shell import single_mie_core_shell
from .material_functions import constant_material, function_material, data_material
from .gmt import spheres, gmt
from .scattering import scattering_per_multipole, absorbption_per_multipole, \
                        extinction_per_multipole, cross_sections