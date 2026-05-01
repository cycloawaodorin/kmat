# Change Log

## 0.0.1 2019/07/19
- First release.

## 0.0.2 2019/12/03
### Bug fix
- Fixed an error in `km_copy_from_work()` for transposed matricies.
Linear algebra methods used with transposed output matricies were affected.
### New feature
- Added` Mat#geo_mean` which returns the geometric mean of self's elements.
- Added some normalizing methods: `Mat#normalize/!`, `Mat#geo_normalize/!` and `Mat#svd_symmetrize/!`.

## 0.0.3 2019/12/03
- Fix a bug in `Mat#geo_normalize`.

## 0.1.0 2022/03/20
- Conformed newer Ruby regulations.
	- Quit using `Random::DEFAULT`.
	- Set `rb_data_type_t::dcompact`.
	- Use `rb_block_call` instead of `rb_iterate`.
	- Use `RB_BLOCK_CALL_FUNC_ARGLIST`.
- Fixed `int`/`size_t` confusion.
- Fixed a bug that `Mat#ge_evd` did not return correct right eigen vectors.
- Fixed some minor bugs.

## 0.1.1 2026/01/07
- Deprecate the global variable `$MatRandom`.
	- It is replaced by a class variable `@@random` under the `Mat` class.
	- The class variable can be obtained by `Mat.random` and can be substituted by `Mat.random=(new_random_instance)`.

## 0.1.2 2026/01/08
- Fixed a bug that `Numeric#*(Mat)` raised TypeError on Ruby 3.4 or later.

## 0.1.3 2026/02/19
- Fixed a bug that there were typos in `Mat#rand_orth`

## 0.1.4 2026/05/01
- Cumulative coding amendment.
