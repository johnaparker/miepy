MiePy
==============
MiePy is a Python module to calculate scattering and absorption properties of spheres and core-shell structures. Mie theory is used, following the procedure in "Absorption and Scattering of Light by Small Particles" by Bohren & Huffman.

The material can be specified by a frequency dependent permittivity and permeability.


Sample Output
--------------
In the image below, scattering intensity is shown for a 100 nm radius dielectric sphere with refractive index n = 3.7. Individual contributions from different multipoles are also shown (eD = electric dipole, etc.).

<p align="center">
  <img src="images/sphere_scattering.png?raw=true" width="600">
</p>


Usage
--------------

For examples and use cases, see examples folder

For full documentation, see docs folder


Installation
--------------
First, clone (or download) this repository and cd into it:

```shell
git clone https://github.com/johnaparker/MiePy && cd MiePy
```

Then, install MiePy via pip (or use just setuptools):

```shell
pip install .                    #using pip
python setup.py install          #using just setuptools
```

To uninstall:

```shell
pip uninstall miepy 
```

To use MiePy without installation, add the MiePy directory to your PYTHONPATH


Dependencies
--------------
For installaltion: setuptools

Required Python modules: numpy, scipy, matplotlib


License
--------------
MiePy is licensed under the terms of the MIT license.
