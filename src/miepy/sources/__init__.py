"""Source module for electromagnetic field sources."""

# Import base classes first to avoid circular imports
# NOTE: Do not reorder - .sources must come before .beams to avoid circular import
from .sources import combined_source as combined_source
from .sources import polarized_propagating_source as polarized_propagating_source
from .sources import propagating_source as propagating_source
from .sources import source as source

# ruff: isort: split
from .beams import beam as beam
from .beams import polarized_beam as polarized_beam
from .beams import reflected_beam as reflected_beam
from .common import azimuthal_beam as azimuthal_beam
from .common import bigaussian_beam as bigaussian_beam
from .common import gaussian_beam as gaussian_beam
from .common import hermite_gaussian_beam as hermite_gaussian_beam
from .common import laguerre_gaussian_beam as laguerre_gaussian_beam
from .common import radial_beam as radial_beam
from .common import shear_beam as shear_beam
from .dft_beams import dft_beam as dft_beam
from .dft_beams import scalar_dft_beam as scalar_dft_beam
from .grid_interpolate import grid_interpolate_source as grid_interpolate_source
from .plane_waves import plane_wave as plane_wave
from .point import point_dipole as point_dipole
from .slm_beam import phase_only_slm as phase_only_slm
from .vsh_sources import vsh_source as vsh_source
