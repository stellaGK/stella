---
title: Input parameters
subtitle: Input parameter namelists for stella
---

# namelist `geo_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`geo_option` | string | `'local'` | Selects the geometry module. Should be one of <ul><li>`local` Miller equilibrium.</li><li>  `vmec` VMEC stellarator equilibrium. Requires VMEC netcdf input file. </li><li> `input.profiles` Reads in General Atomics `input.gacode` file. **This may need to be updated for newer `gacode` files.**</li><li>`miller`same as `local`</li><li>`default` same as `miller`</li></ul> 
`geo_file` | string | `'input.geometry'` | input file used to overwrite selected geometric coefficients below. File uses the same formatting as `.geometry` output.  
`q_as_x` | boolean | `radial_variation` | Uses the safety factor \\( q \\) as the radial coordinate, rather than \\( \psi \\). Used for radially global simulations.
`set_bmag_const` | boolean | `false` | Sets \\( B \\) uniformly to its value at the outboard midplane
`overwrite_bmag` | boolean | `false` | overwrite \\( B \\)
`overwrite_gradpar` | boolean | `false` | overwrite \\( \boldsymbol{b \cdot \nabla b} \\)
`overwrite_gds2` | boolean | `false` |
`overwrite_gds21` | boolean | `false` |
`overwrite_gds22` | boolean | `false` |
 `overwrite_gds23` | boolean | `false` |
 `overwrite_gds24` | boolean | `false` |
 `overwrite_gbdrift` | boolean | `false` |
  `overwrite_gbdrift0` | boolean | `false` |
 `overwrite_cvdrift` | boolean | `false` |
 


# namelist `millergeo_parameters`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`rhoc` | real | 0.5 | radial position of fluxtube \\( r \\)
`rmaj` | real | 3.0 |  device major radius \\( R_0 \\)
`shift` | real | 0.0 | Shafranov shift
`qinp` | real | 1.4 | magnetic safety factor \\( q \\)
`shat` | real | 0.8 | magnetic shear \\( \hat{s} = (r/q) \operatorname{d}q/\operatorname{d}r \\)
`kappa` | real | 0.0 | elongation
`tri` | real | 0.0 | triangularity
`rgeo` | real | 3.0 | geometric center-point of the flux-tube
`betaprim` | real | 0.0 | Normalized pressure gradient \\( \beta' \\)
`kapprim` | real | 0.0 | radial derivative of elongation (for radially global simulation)
`triprim` | real | 0.0 | radial derivative of triangularity (for radially global simulation)
`betadbprim` | real | 0.0 | radial derivative of \\( \beta' \\) (for radially global simulation)
`d2qdr2` | real | 0.0 | second derivative of \\( q' \\) (for radially global simulation)
`d2psidr2` | real | 0.0 | second derivative of \\( \psi \\) (for radially global simulation)
`nzed_local` | integer | 128 | resolution of \\( z \\) on which to calculate geometric coefficients. Splines are then used to interpolate onto the grid used during simulation.
`read_profile_variation` | boolean | `false` | Save information necessary for recalculating Miller equilibrium away from \\( r \\). Used for performing local simulations at different radial locations in order to compare to global simulation.
`write_profile_variation` | boolean | `false` | Recomputes the Miller equilibrium using geometry stored in a file, originally calculated at some \\( r_\mathrm{file} \\), at new \\( r \\).

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
`dkx_over_dky` | real | -1  | Set the desired ratio of \\(\Delta k_x/\Delta k_y\\) for \\(j_\mathrm{twist}\\) when using the stellarator boundary condition.

# namelist `vpamu_grid_parameters`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`nvgrid` | integer | 24 | Number of positive points in \\( v_\parallel \\). Actual grid will be twice as large. Note that \\( v_\parallel = 0\\) is not included.
`nmu` | integer | 12 | Number of points in the magnetic moment \\( \mu \\).
`vpa_max` | real | 3.0 | Maximum \\( v_\parallel \\) on the grid in terms of thermal velocity.
`vperp_max` | real | 3.0 | Maximum \\( v_\perp \\) on the grid in terms of thermal velocity.
`equally_spaced_mu_grid` | boolean | `false` | If false, use Gaussian quadrature points for \\( \mu \\) grid (recommended). Otherwise, use equally spaced grid. Both options exclude \\( \mu = 0 \\).
`conservative_wgts_vpa` | boolean | `false` | Use density-conserving weights for \\( v_\parallel \\). Used primarily by Fokker-Planck collision operator.


# namelist `kt_grids_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`grid_option` | string | `'default'`| Sets the layout of the perpendicular grid. Should be one of <ul><li>`range` use a set range of wavenumbers. Suitable for linear simulation.</li> <li>`box` use a physical box in coordinate space, i.e., a wavenumber grid that satisfies \\(A_\boldsymbol{k} = A^\ast_{-\boldsymbol{k}} \\). *Required* for nonlinear simulation.</li><li>`default` same as `range`</li><li>`annulus` same as `box`</li><li>`nonlinear` same as `box`</li></ul>

# namelist `kt_grids_range_parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`naky` | integer | 1  | Number of points along the binormal mode grid direction.
`nakx` | integer | 1  | Number of points along the radial mode grid direction.
`aky_min` | real | 0.0  | Minimum binormal wavenumber, or binormal wavenumber if `naky` is 1.
`aky_max` | real | 0.0 | Maximum binormal wavenumber if `naky` larger than 1.
`akx_min` | real | 0.0  | Minimum radial wavenumber. Unused if \\( \theta_0 \\) is properly specified.
`akx_max` | real | -1.0 | Maximum radial waveumner.
`theta0_min` | real | 0.0  | Minimum \\( \theta_0 = k_x / \hat{s} k_y \\). Unused if \\( \theta_{0,\textrm{ min}} > \theta_{0,\textrm{ max}} \\).
`theta0_max` | real | -1.0  | Maximum \\( \theta_0 \\).

# namelist `kt_grids_box_parameters`



Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`nx` | integer | 1  | Real space radial resolution.  Resolved number of modes is 2/3 after dealiasing. **Note: checkerboard mode is always excluded.**
`ny` | integer | 1  | Real space binormal resolution.  Resolved number of modes is 2/3 after dealiasing. **Note: checkerboard mode is always excluded.**
`jtwist` | integer | \\( \mathrm{round}(2 \pi \hat{s}) \\)  | Number of independent ballooning modes at \\( k_{y0} \\).
`jtwist_fac` | real | 1.0  | Sets `jtwist` to `jtwist_fac` times the default value if `jtwist` is not specified.
`phase_shift_angle` | real | 0.0  | Binormal phase shift when crossing the parallel boundary. For use with purely periodic boundary conditions, and should have no *statistical* effect for local simulations employing twist-and-shift boundaries. Overwritten for radially global simulations.
`randomize_phase_shift` | boolean | `false`  |  Randomize the binormal phase shift. Unused for radially global simulation.
`x0` | real | `y0`  | Radial extent of the simulation domain. Ignored when using twist-and-shift (linked) boundary conditions.
`y0` | real | -1.0  | Binormal extent of the simulation domain. Overwritten if full-flux-surface is used; **otherwise, must be set.**
`nalpha` | real | 1  | Number of global \\( \alpha \\) grid points. For full-flux-surface simulations.
`centered_in_rho` | boolean | `true` | Ensure domain radially centered around \\( r_0 \\). For use with radially global simulations.
`periodic_variation` | boolean | `false` | Utilize periodic triangle wave radial variation. For use with radially global simulation.

# namelist `layouts_knobs`

These options control the ordering of the coordinate and velocity variables when the domain is decomposed. **At the moment, there is no reason to use anything other than the default options.**

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`xyzs_layout`  | string | `"yxzs"` | Order of the coordinate variables when decomposed for multicore simulations.
`vms_layout` | string | `"vms"` | Order of the velocity variables when decomposed for multicore simulations.


# namelist `physics_flags`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`full_flux_surface` | boolean | `false` | Enables full-flux-surface simulation. Should use finite \\( \rho_\ast \\).
`radial_variation` | boolean | `false` | Enables radially global simulation. Should use finite \\( \rho_\ast \\).
`include_parallel_nonlinearity` | boolean | `false` | Include the parallel nonlinearity. Should use finite \\( \rho_\ast \\).
`include_parallel_streaming` | boolean | `true` | Include parallel streaming.
`include_mirror` | boolean | `true` | Include the mirror term.
`nonlinear` | boolean | `false` | Include the \\( \boldsymbol{E \times B} \\) nonlinear term. **Requires `box` grid.**
`adiabatic_option` | boolean | `false` | Controls how the adiabatic species is treated in quasineutrality. Should be one of <ul><li>`no-field-line-average-term` ion adiabatic species, i.e., \\( \delta n_\mathrm{i}/ n = (Z e/T_\mathrm{i}) \varphi \\). </li> <li> `field-line-average-term` electron adiabatic species, i.e., \\( \delta n_\mathrm{e}/ n = (e/T_\mathrm{e}) (\varphi - \langle \varphi \rangle_\psi) \\).</li><li>`default` same as `no-field-line-average-term`</li><li>`iphi00=0` same as `no-field-line-average-term`</li><li>`iphi00=1` same as `no-field-line-average-term`</li><li>`iphi00=2` same as `field-line-average-term`</li></ul>
`const_alpha_geo` | boolean | `false` | 
`include_pressure_variation` | boolean | `true` | Include kinetic profile variation when running radially global simulations. **Setting this to `false` does *not* currently turn everything off.**
`include_geometric_variation` | boolean | `true` | Include magnetic geomety profile variation when running radially global simulations. **Setting this to `false` does *not* currently turn everything off.**

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


