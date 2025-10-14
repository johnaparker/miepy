MiePy
==============
MiePy is a Python module for the generalized multiparticle Mie theory (GMMT), also known as the aggregate T-matrix method.
MiePy solves the electrodynamics of a collection of spherical or non-spherical scatterers with an arbitrary incident source.

<p align="center">
  <img src="./docs/miepy_thumbnail.png" width="500" alt="Electric field visualization">
  <br>
  <em>Electric field around a 37 particle cluster</em>
</p>

<p align="center">
  <img src="https://jparker.nyc3.digitaloceanspaces.com/gallery/em_three_np_contours.png" width="500" alt="Three particle system">
  <br>
  <em>3D electric field contours around three metal nanoparticles</em>
</p>

Features
--------------
+ **Non-spherical particles** using the T-matrix formulation via the null-field method with discrete sources (NFM-DS). Includes cylinders, spheroids, ellipsoids, cubes and polygonal prisms
+ **Arbitrary incident sources** (plane waves, Gaussian beams, HG and LG beams, point dipoles)
+ Evaluation of cluster **cross-sections** and **optical force and torque** on individual particles
+ **Periodic boundary conditions** with various lattice types (square, hexagonal, etc.) and **mirror and discrete rotational symmetries** for faster calculations
+ Optional **planar interface (substrate)** 
+ **3D scene visualization** using the VPython library
+ Image clusters using a **simulated microscope**
+ **OpenMP parallelization** for systems with larger numbers of particles

Installation
--------------
```shell
pip install miepy
```

If using `uv`:
```
uv venv --python 3.13
uv pip install miepy
source .venv/bin/activate
```

Usage
--------------

See the [examples](src/miepy/examples) folder for how to use MiePy.

Run any of the available [examples](src/miepy/examples) without explicit installation using `uv`:

| Command | Description |
|---------|-------------|
| `uvx miepy dielectric_sphere` | Dielectric sphere scattering and cross-sections |
| `uvx miepy ag_sphere` | Silver sphere scattering and absorption |
| `uvx miepy ag_shell` | Core-shell particle scattering |
| `uvx miepy vary_index` | Scattering intensity vs wavelength and refractive index |
| `uvx miepy fields` | Electric and magnetic field visualization |
| `uvx miepy dimer_scattering` | Au dimer cross-sections |
| `uvx miepy dimer_force` | Force and torque on dimer particles |
| `uvx miepy far_field` | Far-field radiation patterns |
| `uvx miepy whispering_gallery` | Whispering gallery modes in dielectric sphere |
| `uvx miepy focused_gaussian` | Focused Gaussian beam with orbital angular momentum |
| `uvx miepy imaging` | Near-field, far-field, and microscope imaging |

For an overview of the theory, see [docs](./docs) folder.

Install from source
--------------
MiePy uses [vcpkg](https://vcpkg.io/) for C++ dependency management and uv for Python management, which simplifies building across platforms.

**Prerequisites:**
+ GCC and GFORTRAN
+ uv

**Build steps:**

1. Clone MiePy and its submodules:
```shell
git clone https://github.com/johnaparker/miepy.git miepy --recurse-submodules && cd miepy
```

2. Bootstrap vcpkg (first time only):
```shell
./vcpkg/bootstrap-vcpkg.sh
```

3. Install MiePy using uv:
```shell
uv sync
```

4. Optionally, run the tests to verify correctness:
```shell
uv run pytest tests
```

License
--------------
MiePy is licensed under the terms of the GPLv3 license.
