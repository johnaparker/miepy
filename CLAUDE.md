# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands

### Building from Source
The project uses `uv` for package management and `scikit-build-core` for building:

```bash
# Install in development mode
uv pip install -e .

# Build the project (compiles C++ extensions and Fortran executable)
uv build

# Run all tests
pytest tests

# Run a specific test file
pytest tests/test_cluster.py

# Run a specific test function
pytest tests/test_cluster.py::test_off_center_particle
```

### Build System Details
- **Multi-language build**: Python + C++ (via pybind11) + Fortran (NFM-DS T-matrix)
- CMake orchestrates building C++ extensions and Fortran executable
- Material database is automatically extracted from `src/miepy/materials/database.zip` during build
- Fortran executable (`tmatrix`) is built from `src/miepy/tmatrix/nfmds/` and installed to `miepy/bin/`

### Dependencies
Required for building from source:
- CMake (C++ build system)
- Eigen (C++ linear algebra library)
- GNU Scientific Library (GSL)
- GCC and GFORTRAN
- Python 3.11+

## Code Architecture

### Core Solver Classes

MiePy provides two cluster implementations with different optimizations:

1. **`sphere_cluster`** (src/miepy/sphere_cluster.py)
   - Optimized for spherical particles only
   - Uses Mie coefficients directly: `mie_scat[N,2,lmax]` and `mie_int[N,2,lmax]`
   - Faster for sphere-only systems
   - Aggregate T-matrix built via `sphere_aggregate_tmatrix()`

2. **`cluster`** (src/miepy/cluster.py)
   - General implementation for arbitrary particle shapes
   - Uses full T-matrices: `tmatrix[N,2,rmax,2,rmax]`
   - Supports non-spherical particles (spheroids, cylinders, ellipsoids, prisms, etc.)
   - Caches T-matrix calculations for identical particles
   - Aggregate T-matrix built via `particle_aggregate_tmatrix()`

Both classes share a common interface for computing fields, forces, cross-sections, and other observables.

### Vector Spherical Harmonics (VSH) System

The VSH module (src/miepy/vsh/) is the mathematical foundation of the solver:

- **Expansion coefficients**: Fields are expanded in VSH basis with coefficients:
  - `p_src[N,2,rmax]`: Source field expansion at each particle
  - `p_inc[N,2,rmax]`: Incident field (source + scattered from other particles)
  - `p_scat[N,2,rmax]`: Scattered field from each particle
  - `p_int[N,2,rmax]`: Interior field inside each particle
  - `p_cluster[2,rmax]`: Cluster scattering coefficients around system origin

- **Index conventions**:
  - `lmax`: Maximum angular momentum order
  - `rmax = lmax*(lmax + 2)`: Total number of modes per particle/polarization
  - First dimension [2]: Electric (0) and magnetic (1) multipoles
  - Second dimension [rmax]: Combined (n,m) mode index via `mode_indices(lmax)`

- **VSH modes**: Defined by `vsh_mode` enum (incident, outgoing, interior)

### Interaction Matrix System

For multi-particle systems with `interactions=True`:

1. **Source decomposition** (`_solve_source_decomposition`):
   - Expands incident source in VSH basis at each particle position
   - Stores in `p_src`

2. **Linear system** (`_solve_interactions`):
   - Builds aggregate T-matrix encoding particle scattering + inter-particle coupling
   - Solves: `p_inc = (I + T_agg)^(-1) * p_src` using BiCGSTAB iterative solver
   - T-matrix assembly happens in C++ for performance

3. **Scattering coefficients**:
   - For spheres: `p_scat = p_inc * mie_scat` (element-wise)
   - For general particles: `p_scat = einsum('naibj,nbj->nai', tmatrix, p_inc)`

### Particle System

All particles inherit from `particle_base` (src/miepy/particles/particle_base.py):

- **Spherical particles**: Use analytical Mie theory (fast)
- **Non-spherical particles**: Compute T-matrix via NFM-DS Fortran code
  - Creates input files in temporary directory
  - Calls external `tmatrix` executable
  - Parses output files for T-matrix data
  - Caching via `_dict_key()` prevents redundant calculations

Available particle types (src/miepy/particles/):
- `sphere`, `core_shell`: Analytical Mie solution
- `spheroid`, `cylinder`, `ellipsoid`: Axisymmetric NFM-DS
- `regular_prism`, `cube`: Non-axisymmetric NFM-DS
- `sphere_cluster_particle`: Cluster treated as single particle

### Source System

Sources (src/miepy/sources/) define incident electromagnetic fields:

- **Plane waves**: `plane_wave.from_string(polarization='x')` or custom k-vector
- **Gaussian beams**: Paraxial and focused implementations
- **Point sources**: Point dipoles for LDOS calculations
- **Structured beams**: Hermite-Gaussian, Laguerre-Gaussian via DFT decomposition

All sources implement:
- `E_field(x, y, z, k)` and `H_field(x, y, z, k)`: Direct field evaluation
- `structure(positions, k, lmax)`: VSH decomposition at particle positions

### Material System

Materials (src/miepy/materials/ and src/miepy/material_functions/):

- **Predefined materials**: `miepy.materials.Ag()`, `miepy.materials.Au()`, etc.
- **Database**: Refractive index data in `materials/database/` (extracted from zip)
- **Custom materials**:
  - `constant_material(eps, mu)`: Frequency-independent
  - `function_material(eps_func, mu_func)`: Functional form
  - `data_material(wavelengths, n_data, k_data)`: Interpolated data

### C++ Backend

Performance-critical operations in cpp/src/:

- VSH translations (particle-particle coupling)
- Aggregate T-matrix assembly
- Linear system solving
- Special function evaluations

Exposed via pybind11 as `miepy.cpp` module.

## Common Development Patterns

### Computing Fields

```python
# Total field (source + scattered)
E = cluster.E_field(x, y, z)

# Scattered field only
E_scat = cluster.E_field(x, y, z, source=False)

# Far-field
E_far = cluster.E_field(r, theta, phi, far=True, spherical=True)
```

### Cross-sections

```python
# Total cluster cross-sections
C = cluster.cross_sections()  # Returns named tuple (scattering, absorption, extinction)

# Per-particle cross-sections
C_i = cluster.cross_sections_of_particle(i)

# Per-multipole (for cluster around origin)
C_multipole = cluster.cross_sections_per_multipole()
```

### Updating Geometry

```python
# For sphere_cluster
cluster.update_position(new_positions)

# For cluster (supports orientation changes)
cluster.update(position=new_positions, orientation=new_orientations)
```

## Important Implementation Details

### Coordinate Systems
- Default: Cartesian (x, y, z)
- Spherical: (r, theta, phi) with theta from z-axis, phi from x-axis
- Many methods support `spherical=True` parameter for I/O in spherical coordinates
- VSH expansions are naturally in spherical coordinates

### Origins
- Cluster has a system `origin` (default [0,0,0])
- Use `origin='auto'` to automatically center on particle geometry
- Cross-sections are origin-independent, but cluster coefficients (`p_cluster`) depend on origin choice

### Symmetries and Interfaces
- Periodic symmetries (lattices) supported in `sphere_cluster` only
- Interface (planar substrate) partially implemented but marked as not fully supported in `cluster`

### Performance Considerations
- T-matrix calculations for non-spherical particles are expensive (Fortran executable call)
- Identical particles share T-matrix via caching (based on `_dict_key()`)
- C++ backend used for interaction matrix assembly and field evaluations
- For large systems, BiCGSTAB convergence can be slow; consider reducing `lmax`
