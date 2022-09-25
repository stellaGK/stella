---
title: Input parameters
subtitle: Input parameter namelists for stella
---

# namelist `geo_knobs`
# namelist `millergeo_parameters`
# namelist `vmec_parameters`

# namelist `zgrid_parameters`

Used to control the resolution and boundary condition of the parallel coordinate.


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`nzed` | integer | 24 | description
`nperiod` | integer | 1 | Controls number of poloidal turns in the flux tube, which is given by \\(n_\mathrm{pol} = 2 n_\mathrm{period} -1\\)
`ntubes` | integer | 1 | Number of cars in the flux tube train
`boundary_option` | string | `'default'` | Parallel boundary condition, should be one of <ul><li>`default` zero incoming, no radial connections, suitable for linear simuatlions</li><li>  `linked` zero incoming with twist-and-shift boundary condition, suitable for nonlinear simulations</li><li>  `stellarator` Martin's twist-and-shift boundary for stellarators with local shear</li><li> `periodic` periodic for all modes, with option of a phase shift across the parallel boundary</li><li>`zero`same as `default`</li><li>`unconnected` same as `default`</li><li>`self-periodic` same as `periodic`</li></ul> Note that for all of these options, zonal modes are treated as periodic.
`zed_equal_arc` | boolean | `false` | Arc length is used as the parallel coordinate, so that  \\( \boldsymbol{\hat{b} \cdot \nabla} \theta \\) is constant.
`shat_zero` | real | \\(10^{-5}\\) | The minimum value of magnetic shear (\\(\hat{s}\\)) needed to employ the linked boundary condition. Otherwise, periodic boundary conditions are used
`grad_x_grad_y_zero` | real | \\(10^{-5}\\) | minimum \\(\nabla x \cdot \nabla y \\) value at the end of the flux-tube which we assume periodic boundary conditions instead of the stellarator symmetric ones when using stellarator boundary condition
`dkx_over_dky` | real | -1  | Set the desired ratio of \\(dkx/dky\\) for \\(j_\mathrm{twist}\\) when using the stellarator boundary condition.

# namelist `vpamu_grid_parameters`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
*variable* | integer | 1 | description


# namelist `kt_grids_knobs`
# namelist `kt_grids_range_parameters`
# namelist `kt_grids_box_parameters`
# namelist `layouts_knobs`

# namelist `physics_flags`
# namelist `parameters`
# namelist `knobs`
# namelist `init_g_knobs`
# namelist `species_knobs`
# namelist `species_parameters`
# namelist `time_advance_knobs`
# namelist `stella_diagnostics_knobs`

# namelist `multibox_parameters`
# namelist `sources`


# namelist `dissipation`
# namelist `collisions_dougherty`
# namelist `collisions_fp`
# namelist `hyper`

# namelist `sfincs_input`
# namelist `neoclassical_input`


