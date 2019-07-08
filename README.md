MiePy
==============
MiePy is a Python module for the generalized multiparticle Mie theory (GMMT), also known as the aggregate T-matrix method. MiePy solves the electrodynamics of a collection of spherical or non-spherical scatterers with an arbitrary incident source.

Features
--------------
+ **Non-spherical particles** using the T-matrix formulation via the null-field method with discrete sources (NFM-DS). Includes cylinders, spheroids, ellipsoids, cubes and polygonal prisms
+ **Arbitrary incident sources** (plane waves, Gaussian beams, HG and LG beams, point dipoles)
+ Evaluation of cluster **cross-sections** and **optical force and torque** on individual particles
+ **Periodic boundary conditions** with various lattice types (square, hexagonal, etc.) and **mirror and discrete rotational symmetries** for faster calculations
+ Optional **planar interface (substrate)** 
+ **3D scene visualization** using the the VPython library
+ Image clusters using a **simulated microscope**
+ **OpenMP parallelization** for systems with larger numbers of particles

Usage
--------------

For examples and use cases, see examples folder.

For an overview of the theory, see docs folder.

Installation
--------------
If NumPy is not already installed, it must be installed prior to MiePy's installation
```shell
pip install numpy
```
Then install MiePy
```shell
pip install miepy
```
MiePy is also available via Conda
```shell
conda install -c japarker miepy
```

Install from source
--------------
To build MiePy from source, first install the required dependencies:

+ [CMake](https://cmake.org/install/) (C++ build system)
+ [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (C++ linear algebra library)
+ [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
+ GCC and GFORTRAN
+ Python 3 and pip

Then, install MiePy using pip
```shell
pip install miepy --no-binary
```

To build the latest development version, clone MiePy and its submodules:
```shell
git clone https://github.com/johnaparker/miepy.git miepy --recurse-submodules && cd miepy
```
and install MiePy using pip
```shell
pip install .
```
Optionally, run the tests to verify correctness:
```shell
pytest tests
```

License
--------------
MiePy is licensed under the terms of the GPLv3 license.
