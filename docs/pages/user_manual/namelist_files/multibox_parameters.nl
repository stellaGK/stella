# namelist `multibox_parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`boundary_size` | integer | 4 | Number of collocation points in the boundary region.
`krook_size` | integer | 0 | Number of collocation points in the Krook region. Will automatically max out at `boundary_size`.
`zf_option` | string | `'default'` | *Experimental* - Set how the zonal mode is handled during the communication of the multiple-flux-tube boundary condition. Should be one of
`krook_option` | string | `''` | Shape of the Krook operator in the Krook region.
`RK_step` | boolean | `false` | Communicate the multiple-flux-tube boundary condition at every implicit and Runge-Kutta substep, rather than one per timestep. 
`nu_krook_mb` | real | 0.0 | Strength (i.e. damping rate) of the Krook operator in the boundary region.
`mb_debug_step` | integer | 1000 | 
`krook_exponent` | real | 0.0 |  
`comm_at_init` | boolean | `false` | 
`phi_bound` | integer | 0 | 
`phi_pow` | real | 0.0 | 
`krook_efold` | real | 3.0 | 
`use_dirichlet_BC` | boolean | `false` | 
`LR_debug_option` | string | `'default'`| 
`smooth_ZFs` | boolean | `false` | *Experimental* - Try smoothing the zonal flows across the interface between boundary and physical region. *Never really worked properly.* 
