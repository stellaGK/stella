# namelist `vmec_parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`alpha0` | real | 0.0 |
`zeta_center` | real | 0.0  |
`nfields_periods` | real | -1.0 |
`torflux` | real | 0.6354167|
`zgrid_scalefac` | real | 1.0 |
`zgrid_refinement_factor` | real | 4 if `zed_equal_arc`, 1 otherwise|
`surface_option` | integer | 0 |
`verbose` | boolean | `true` |
`vmec_filename` | string | `'equilibria/wout_w7x_standardConfig.nc'` |
`gradpar_zeta_prefac` | real | 1.0 |
`n_tolerated_test_arrays_inconsistencies` | integer | 0 | Non-negative integer. Ignores a this number of test_arrays inconsistencies. Useful when running with low-resolution VMEC files, which can be useful for optimization. 0 or 1 are recommended values.
