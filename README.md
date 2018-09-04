MiePy
==============
MiePy is a Python module for the generalized multiparticle Mie theory (GMMT), also known as the aggregate T-matrix method. MiePy solves the electrodynamics of a collection of spherical or non-spherical scatterers with an arbitrary incident source.

Features
--------------
+ **Non-spherical particles** using the T-matrix formulation via the null-field method with discrete sources (NFM-DS). Includes **cylinders, spheroids, ellipsoids, cubes and polygonal prisms
+ **Arbitrary incident sources** (plane waves, Gaussian beams, HG and LG beams, point dipoles) via near-field point matching or far-field integration
+ **Periodic boundary conditions** with various lattice types (square, hexagonal, etc.) and **mirror and discrete rotational symmetries** for faster calculations
+ Evaluation of cluster **cross-sections** and **optical force and torque** on individual particles
+ **OpenMP parallelization** for systems with larger numbers of particles

Usage
--------------

For examples and use cases, see examples folder.

For full documentation, see docs folder.

Installation
--------------
Install required dependencies:

+ [CMake](https://cmake.org/install/) (C++ build system)
+ [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (C++ linear algebra library)
+ [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
+ GCC and GFORTRAN
+ Python 3 and pip

Clone MiePy and its submodules:
```shell
git clone https://johnapark@bitbucket.org/johnapark/miepy.git --recurse-submodules && cd miepy
```

Build and install MiePy using pip:
```shell
python setup.py build_ext
pip install . -e --user
```

Optionally, run tests to verify correctness:
```shell
pytest tests
```

Sample Output
--------------

License
--------------
MiePy is licensed under the terms of the MIT license.
