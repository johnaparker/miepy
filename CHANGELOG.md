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
