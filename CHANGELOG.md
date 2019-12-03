# Change Log

## 0.0.1 2019/07/19
- First release.

## 0.0.2 2019/12/03
### Bug fix
- Fixed an error in km_copy_from_work() for transposed matricies.
Linear algebra methods used with transposed output matricies were affected.
### New feature
- Added Mat#geo_mean which returns the geometric mean of self's elements.
- Added some normalizing methods: Mat#normalize/!, Mat#geo_normalize/! and Mat#svd_symmetrize/!
