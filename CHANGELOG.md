## 0.4 - 2019-05-30
- Support for including an interface (substrate) in the simulation; all sources can be transmitted and reflected
- Non-axisymmetric particles are now supported; ellipsoids, cubes, and regular prisms are defined
- Focal fields of focused beams are computed exactly using the angular spectrum rather than via paraxial limit expressions
- 3D scene visualization using the VPython library
- Microscope module for producing simulated images of clusters through a lens system
- New metal material for perfectly conducting objects
- Method to compute the local density of states (LDOS)
- Faster T-matrix initialization for particles of the same geometry by caching the non-rotated T-matrix
- Performance critical code re-written in C++ and binded using pybind11; large speed-up for smaller clusters
- Experimental support for OpenMP for shared-memory parallel execution
- Platform support: now builds on Linux, MacOS and Windows

## 0.3 - 2018-07-16
- Extension to non-spherical particles (spheroids and cylinders) using the T-matrix method
- Plane waves and beams now have theta and phi variables to control their incident direction
- 30x performance improvement for larger clusters

## 0.2 - 2018-06-13
- New cluster class and removal of old spheres class for easier cluster creation
- More useful E and H field functions that include interior field expansion, an option to mask the interior fields, a far-field option for more efficient far-field expansion, and a spherical option for evaluation the fields in spherical coordinates
- Paraxial and non-paraxial focused beams have been added using a new source decomposition method (a hybrid near-field point-matching and far-field integration method)
- All beams can now be power and phase controlled
- New beams: bigaussian and paraxial (a generic paraxial beam with an arbitrary scalar potential)
- VSH expansion coefficients combined into single array instead of two throughout the API
- New example: focused Gaussian beam
- Removed old cross-sections module

## 0.1 - 2018-03-29
