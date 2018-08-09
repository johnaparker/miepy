from . import special
from . import old_special

from .mode_indices import rmax_to_lmax, lmax_to_rmax, mode_indices
from .vsh_functions import (Emn, vsh_mode, get_zn, VSH, vsh_normalization_values,
                             get_zn_far, VSH_far, vsh_normalization_values_far)
from .vsh_translation import vsh_translation 
from .vsh_rotation import vsh_rotation_matrix, rotate_expansion_coefficients
from .expansion import expand_E, expand_E_far, expand_H, expand_H_far
from .decomposition import (near_field_point_matching, far_field_point_matching, 
                            integral_project_fields_onto, integral_project_fields,
                            integral_project_source, integral_project_source)
from .cluster_coefficients import cluster_coefficients
