from . import old_special, special

# Import these last to avoid circular import issues
# (they import from miepy.vsh and use symbols defined above)
from .cluster_coefficients import cluster_coefficients
from .decomposition import (
    far_field_point_matching,
    integral_project_fields,
    integral_project_fields_onto,
    integral_project_source,
    near_field_point_matching,
)
from .expansion import expand_E, expand_E_far, expand_H, expand_H_far
from .mode_indices import lmax_to_rmax, mode_indices, rmax_to_lmax
from .vsh_functions import (
    VSH,
    Emn,
    VSH_far,
    get_zn,
    get_zn_far,
    vsh_mode,
    vsh_normalization_values,
    vsh_normalization_values_far,
)
from .vsh_rotation import rotate_expansion_coefficients, vsh_rotation_matrix
from .vsh_translation import vsh_translation
