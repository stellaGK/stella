---
title: Namelist parameters
subtitle: Input parameter namelists for stella
---

# namelist `geo_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`geo_option` | string | `'local'` | Selects the geometry module. Should be one of <ul><li>`local` Miller equilibrium.</li><li>  `vmec` VMEC stellarator equilibrium. Requires VMEC netcdf input file. </li><li> `input.profiles` Reads in General Atomics `input.gacode` file. **This may need to be updated for newer `gacode` files.**</li><li>`miller`same as `local`</li><li>`default` same as `miller`</li></ul> 
`geo_file` | string | `'input.geometry'` | input file used to overwrite selected geometric coefficients below. File uses the same formatting as `.geometry` output.  
`q_as_x` | boolean | `radial_variation` | Uses the safety factor \\( q \\) as the radial coordinate, rather than \\( \psi \\). Used for radially global simulations.
`set_bmag_const` | boolean | `false` | Sets \\( B \\) uniformly to its value at the outboard midplane
`overwrite_bmag` | boolean | `false` | overwrite \\( B \\).
`overwrite_gradpar` | boolean | `false` | overwrite \\( \boldsymbol{ b \cdot \nabla z} \\).
`overwrite_gds2` | boolean | `false` | overwrite \\(  \lvert \nabla \alpha \rvert^2 (\textrm{d}\psi / \textrm{d}r)^2 \\).
`overwrite_gds21` | boolean | `false` | overwrite \\( \lvert \nabla q \cdot \nabla \alpha \rvert (\textrm{d}\psi / \textrm{d}r)^2 \\).
`overwrite_gds22` | boolean | `false` | overwrite \\( \lvert \nabla q \rvert^2 (\textrm{d}\psi / \textrm{d}r)^2 \\).
 `overwrite_gds23` | boolean | `false` | overwrite \\(  \nabla \theta \cdot [\nabla \alpha \times (\nabla r \times \nabla \alpha)] (\textrm{d}\psi / \textrm{d}r)^2 / B^2 \\).
 `overwrite_gds24` | boolean | `false` | overwrite \\( \nabla \theta \cdot [\nabla r \times (\nabla r \times \nabla \alpha)] (\textrm{d}\psi / \textrm{d}r)^2 (q/r) / B^2\\).
 `overwrite_gbdrift` | boolean | `false` | overwrite \\( 2 (\boldsymbol{b} \times \nabla B \cdot \nabla \alpha) (\textrm{d}\psi / \textrm{d}r) /B^2 \\).
  `overwrite_gbdrift0` | boolean | `false` | overwrite \\( 2 (\boldsymbol{b} \times \nabla B \cdot \nabla q) (\textrm{d}\psi / \textrm{d}r) /B^2  \\).
 `overwrite_cvdrift` | boolean | `false` | overwrite \\( 2(\textrm{d}\psi / \textrm{d}r)[\boldsymbol{b} \times (\boldsymbol{b \cdot \nabla b})]\cdot \nabla \alpha / B \\)


# namelist `millergeo_parameters`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`rhoc` | real | 0.5 | radial position of fluxtube \\( r \\)
`rmaj` | real | 3.0 |  device major radius \\( R_0 \\)
`shift` | real | 0.0 | Shafranov shift
`qinp` | real | 1.4 | magnetic safety factor \\( q \\)
`shat` | real | 0.8 | magnetic shear \\( \hat{s} = (r/q) \textrm{d}q/\textrm{d}r \\)
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
`n_tolerated_test_arrays_inconsistencies` | integer | 0 | Non-negative integer. Ignores a this number of test_arrays inconsistencies. Useful when running with low-resolution VMEC files, which can be useful for optimization. 0 or 1 are recommended values.


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
`kyspacing_option` | string | `'default'`  | Sets the spacing between \\(a k_y\\) grid points, available options are : <ul><li>`default` same as `linear`</li><li>`linear` linear spacing in \\(a k_y\\)</li><li>`exponential` linear spacing in \\(\log{(a k_y)}\\)</ul>


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

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`beta`  | float | 0.0 | Plasma \\( \beta \\). Currently has no effect.
`zeff`  | float | 1.0 | Effective charge number for use with *effective* electron-ion and electron-impurity collisions in the Fokker-Planck collision operator (see `ecoll_zeff`).
`tite`  | float | 1.0 | Ratio of ion to electron temperature, \\( T_\mathrm{i}/T_\mathrm{e} \\). Used in quasineutrality when adiabatic species is used.
`nine`  | float |1.0 | Ratio of ion to electron density, \\( n_\mathrm{i}/n_\mathrm{e} \\). Used in quasineutrality when adiabatic species is used.
`rhostar`  | real | -1.0 | The gyrokinetic expansion parameter \\( \rho_\mathrm{th,ref}/a_\mathrm{ref} \\). For effects beyond the flux-tube limit (full-flux-surface, radially global, neoclassical terms, etc...). Overwritten if `irhostar` is positive.
`irhostar`  | real | -1.0 | Sets `rhostar = 1.0 / irhostar` if positive.
`vnew_ref`  | real | -1.0 | Reference collision frequency. Various input options will overwrite this if it is negative.
`g_exb`  | real | 1.0 | Equilibrium \\( \boldsymbol{E \times B} \\) shear rate. More specifically, \\( \gamma_\boldsymbol{ E \times B} = (r/q) (\textrm{d}\omega / \textrm{d}r) R_0/\sqrt{2}v_\mathrm{th,ref}\\). Uses the Hammett wavenumber shift method, with nonlinear corrections proposed by McMillan.
`g_exbfac`  | real | 1.0 | Prefactor for perpendicular component of equilibrium \\( \boldsymbol{E \times B} \\) flow shear. Setting to 0.0 turns this component off.
`omprimfac`  | real | 1.0 | Prefactor for parallel component of equilibrium \\( \boldsymbol{E \times B} \\) flow shear. Setting to 0.0 turns this component off.


# namelist `knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`fphi`  | real | 1.0 | Prefactor for electrostatic potential \\( \varphi \\) wherever it appears. 
`fapar`  | real | 1.0 | Prefactor for fluctuating vector potential \\( A_\parallel \\) wherever it appears. *Currently has no effect*.
`fbpar`  | real | -1.0 | Prefactor for fluctuating parallel magnetic field \\( B_\parallel \\) wherever it appears. *Currently has no effect*.
`delt`  | real | 0.0 | Initial simulation timestep. CFL constraints may change this throughout the simulation.
`nstep`  | integer | -1 | Number of simulation timesteps.
`tend`  | real | -1.0 | End-time of the simulation. If not set, then not used.
`delt_option`  | string | '`check_restart `' | How to handle setting the timestep on restart. Should be one of <ul><li>`check_restart` automatically checks the restart file for last saved time step.</li><li>`set_by_hand` use `delt` from input file.</li><li>`default` same as `check_restart`.</li></ul>
`lu_option`  | string | `default` | Parallelization of the LU decomposition. Should be one of <ul><li>`none` no parallelization.</li><li>`none` same as `default`</li><li>`local` Parallelized locally on a core using shared memory. Best case speed-up is \\(j_\mathrm{twist}\times (\textrm{cores per node}). \\).</li><li>`global` parallelized over all cores. Currently only works on experimental branch `development/pLU_scalapack` which uses ScaLAPACK. </li></ul>If compiled with `HAS_ISO_C_BINDING`, then it is ***strongly*** recommended to run with `lu_option='local'`.
`avail_cpu_time`  |real |  \\( 10^{10} \\) | Available CPU time **in seconds**. Useful for cleanly ending a run before allocated time runs out.
`cfl_cushion_upper`  | real | 0.5 | Safety margin for the CFL condition. 
`cfl_cushion_middle`  | real | 0.25 | Safety margin used for setting a new *smaller* timestep when the timestep exceeds cfl_cushion_upper \* CLF_dt, or for setting a new "larger" timestep when the timesetps is smaller than cfl_cushion_lower \* CFL_dt. 
`cfl_cushion_lower`  | real |  0.00001 | Lowest time-step based on CFL condition. 
`delt_max`  | real | -1 | If positive, then set the maximum timestep to `delt_max`; otherwise, the maximum time step will be the initial one.
`stream_implicit`  | boolean | `true` | Calculate parallel streaming implicitly using the response matrix approach.
`mirror_implicit`  | boolean | `true` | Calculate the mirror term implicitly.
`driftkinetic_implicit`  | boolean | `false` | When calculating parallel streaming, only include the non-gyroaveraged electrostatic potential \\( \varphi \\) in the implicit calculation, and calculate the portion resulting from \\( \varphi - \langle \varphi \rangle_\boldsymbol{R} \\) explicitly.
`drifts_implicit`  | boolean | `false` | Calculate the magnetic and \\( \omega_\ast \\) drifts implicitly as an extra term to the operator splitting.
`stream_matrix_inversion`  | boolean | `false` | Use a different tri-diagional solver for parallel streaming. 
`maxwellian_inside_zed_derivative`  | boolean | `false` | *Experimental* - Evaluate the parallel streaming term with the Maxwellian background inside the parallel derivative, and also include the extra term proportional to \\( \partial_z B\\) that results from the product rule.
`mirror_semi_lagrange`  | boolean | `true` | Use semi-Lagrange solve for mirror term. Otherwise, use tri-diagonal matrix solve.
`mirror_linear_interp`  | boolean | `false` | Use linear, rather than fourth-order, interpolation when using semi-Lagrange approach for mirror term.
`zed_upwind`  | real |  0.02 | Amount of spatial upwinding in \\( z \\) when implicit solve is used. Recommended values: 0.02–0.05.
`vpa_upwind`  |  real | 0.02 | Amount of upwinding in \\( v_\parallel \\) when implicit mirror term is used with the matrix solve. Recommended values: 0.02–0.05. *No effect when semi-Lagrange is used.*
`time_upwind`  | real  | 0.02 | Amount of temporal in \\( t \\) when implicit solve is used. Recommended values: 0.02–0.05.
`fields_kxkyz`  | boolean | `false` | Calculate electromagnetic fields with a local velocity grid. **Requires MPI all-to-all redistrution, and so *not* recommended.**
`mat_gen`  | boolean | `false` | Write out response matrices. **`lu_option='local'` makes this obsolete in most cases**.
`mat_read`  | boolean | `false` | Read in response matrices. **`lu_option='local'` makes this obsolete in most cases**.
`rng_seed`  | integer | -1 | Seeds the random number generator used for the `noise` initial condition. If negative, then a seed is generated from the current time.
`ky_solve_radial`  |  integer |  0 | How many \\( k_y \\) modes, starting from the zonal mode, for which quasineutrality is to be calculated exactly, rather than perturbatively, for radially global simulation. **It is recommended to calculate quasineutrality for *all* modes exactly, so set this to naky when running radially global.**
`ky_solve_real`  | boolean | `false` | *Experimental* - solve quasineutrality in *real space excluding the boundary region* when performing a radially global simulation.


# namelist `init_g_knobs`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`ginit_option` | string | `'default'` | Distribution function initialization option. Should be one of <ul><li>`default` Gaussian pulse along the zed dimension for every wavenumber. Suitable for linear simulation. </li><li>`noise` white noise for every wavenumber. Suitable for nonlinear simulation. </li><li>`many` Read from (multiple) restart netcdf files. </li><li>`nltest` Not implemented. </li><li>`kxtest` Not implemented. </li><li>`kpar` Distribution function with particular moments and parallel shapes set by options below. </li><li>`rh` Rosenbluth-Hinton initialization of \\( k_y = 0\\) zonal mode. Useful for testing linear parallel streaming and magnetic drift physics. </li><li>`remap`  Initalizates one or two single eigenmodes with finite amplitude. Used for testing wavenumber remapping implementation of equilibrium \\( \boldsymbol{E \times B} \\) flow shear and, in particular, its nonlinear corrections. </li></ul>
`width0` | real | -3.5 | If using `ginit_option='default'`, sets the width of the Gaussian pulse.
`phiinit` | real | 1.0 | Rough size of the resulting root-mean-square electrostatic potential \\( \varphi \\)from initialization. Actual size will vary depending on species options, `ginit_option`, etc... See also `scale_to_phiinit`.
`scale_to_phiinit` | boolean | `false` | If true, then upon initialization (not restart!) rescale the distribution function so that the root-mean-squared electrostatic potential \\( \varphi \\) is *precisely* `phiinit`, regardless of any other option.
`restart_file` | string | `RUN_NAME.nc` | Name of the restart file. Note that  `RUN_NAME` derives from `RUN_NAME.input`, the input file used to run stella. 
`restart_dir` | string | `./` | Restart directory. *Tip to users:* using `./nc/` is often very convenient.
`read_many` | boolean | `true` | If false, then use MPI_IO to write to one restart file. Otherwise, every core writes to a separate restart file (*recommended*).
`chop_side` | boolean | `false` | Zero out the distribution function on one side of the \\(z\\) grid. See `left`.
`left` | boolean | `true` | If `chop_side` is true, then zero out the distribution on the left size (\\( z < 0\\)) of the zed domain. Otherwise, chop it on the right side (\\( z > 0 \\)).
`scale` | real | 1.0 | Rescale the distribution function by `scale` upon restart.
`tstart` | real | 0.0 | Start time of the simulation. Overwritten when restarting.
`zf_init` | real | 1.0 | Prefactor for the \\( k_y = 0 \\) zonal modes. Setting to 0.0 zeros them out upon initialization.
`den0` | real | 1.0 | Constant density perturbation component when using `ginit_option='kpar'`.
`upar0` | real | 0.0 | Constant parallel velocity perturbation component when using `ginit_option='kpar'`.
`tpar0` | real | 0.0 | Constant parallel temperature perturbation component when using `ginit_option='kpar'`.
`tperp0` |real | 0.0 | Constant perpendicular temperature perturbation component when using `ginit_option='kpar'`.
`den1` | real | 0.0 | \\( \cos(\theta) \\) density perturbation component when using `ginit_option='kpar'`.
`upar1` | real | 0.0 | \\( \cos(\theta) \\) parallel velocity perturbation component when using `ginit_option='kpar'`.
`tpar1` | real | 0.0 | \\( \cos(\theta) \\) parallel temperature perturbation component when using `ginit_option='kpar'`.
`tperp1` | real | 0.0 | \\( \cos(\theta) \\) perpendicular component perturbation component when using `ginit_option='kpar'`.
`den2` | real | 0.0 | \\( \cos(2\theta) \\) density perturbation component when using `ginit_option='kpar'`.
`upar2` | real | 0.0 |  \\( \cos(2\theta) \\) parallel velocity perturbation component when using `ginit_option='kpar'`.
`tpar2` | real | 0.0 |  \\( \cos(2\theta) \\) parallel temperature perturbation component when using `ginit_option='kpar'`.
`tperp2` | real | 0.0 |  \\( \cos(2\theta) \\) perpendicular temperature perturbation component when using `ginit_option='kpar'`.
`refac` | real | 1.0 | Controls the real part of the distribution function when using `ginit_option='kpar'`.
`imfac` | real | 0.0 | Controls the imaginary part of the distribution function when using `ginit_option='kpar'`.
`kxmax` | real | \\( 10^{100} \\)| The maximum radial wavenumber up to which modes are initialized when using `ginit_option='rh'`.
`kxmin` | real | 0.0 | The maximum minimum wavenumber up fro, which modes are initialized when using `ginit_option='rh'`.


# namelist `species_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`nspec` | integer | 2  | Number of species.
`species_option` | string | `'stella'` | How to read in species data. Should be one of <ul><li>`stella` Read from stella input file.</li><li>  `euterpe` Read form euterpe file. </li><li> `input.profiles` Reads in General Atomics `input.gacode` file. **This may need to be updated for newer `gacode` files.**</li><li>`default` same as `stella`</li></ul> 
`read_profile_variation` | boolean | false | Save information necessary for recalculating kinetic profile quantities away from \\( r \\). Used for performing local simulations at different radial locations in order to compare to global simulation.
`write_profile_variation` | boolean | false | Recomputes the kinetic profile information using kinetic profile quantities stored in a file, originally calculated at some \\( r_\mathrm{file} \\), at new \\( r \\).
`ecoll_zeff` | boolean | false | If true, use an effective intra-species electron-ion collision rate using `zeff`.


# namelist `species_parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`z` | integer  | 1 | charge number.
`mass` | real | 1.0 | particle mass, normalized to \\( m_\mathrm{ref} \\). 
`dens` | real | 1.0 | species density, normalized to \\( n_\mathrm{ref} \\). 
`temp` | real | 1.0 | species temperature, normalized to \\( T_\mathrm{ref} \\). 
`fprim` | real | -999.9 | Density gradient scale length \\(a_\mathrm{ref}L_{n_s}^{-1} = - a_\mathrm{ref} \textrm{d} \ln n_s/ \textrm{d}r \\)). Note the negative sign.
`tprim` | real | -999.9 | Temperature gradient scale length \\(a_\mathrm{ref}L_{T_s}^{-1} = - a_\mathrm{ref} \textrm{d} \ln T_s/ \textrm{d}r \\)). Note the negative sign.
`d2ndr2` | real | 0.0 | Second derivative of density, normalized to \\( a_\mathrm{ref}^2/n_\mathrm{ref}\\). For use with radially global simulation.
`d2Tdr2` | real | 0.0 | Second derivative of temperature, normalized to \\( a_\mathrm{ref}^2/T_\mathrm{ref}\\). For use with radially global simulation.
`bess_fac` | real | 1.0 | Perfactor for Bessel function argument. Setting to 0.0 renders particle drift-kinetic.
`type` | string | `'default'` | Particle type. Should be one of <ul><li>`ion` ion species.</li><li>`default` same as `ion`. </li><li>  `electron` electron species. </li><li> `e` same as `electron` </li><li>`beam` slowing down species. </li><li>`fast` same as `beam`. </li><li>`alpha` same as `beam`. </li><li>`slowing-down` same as `beam`. </li><li>`trace` tracer species. </li></ul> 


# namelist `time_advance_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`xdriftknob` | real | 1.0 | Prefactor for radial magnetic drift. Setting to 0.0 turns the term off.
`ydriftknob` | real | 1.0 | Prefactor for the binormal magnetic drift. Setting to 0.0 turns the term off.
`wstarknob` | real | 1.0 | Prefactor for the \\( \omega_\ast \\) term. Setting to 0.0 turns the term off.
`explicit_option` | string | `'default'` | Chooses the Runge-Kutta scheme for the explicit integration. Should be one of <ul><li>`rk2` second-order Runge-Kutta. </li> <li>`rk3` third-order strong-stability-preserving Runge-Kutta (recommended). </li><li>`rk4` fourth-order Runge-Kutta. </li><li>`default` same as `rk3`. </li></ul> Note that higher-order Runge-Kutta schemes can increase memory usage.
`flip_flop` | boolean | `false` | Utilize the flip-flopping approach that flips the integration order every time-step. Should increase time accuracy, at least linearly. Does sometimes lead to spurious oscillations.


# namelist `stella_diagnostics_knobs`

These options control which diagnostics are output by stella. **NOTE: stella safely appends to both ASCII and netCDF files upon restart, so there is no need to copy files and move them around.** Also note that radial fluxes are always written to an ASCII file named `RUN_NAME.fluxes`.

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`nwrite` | integer | 50 | Output cadence (i.e. number of iterations) between output diagnostics.
`navg` | integer | 50 | Number of timesteps over which to average the real and imaginary mode frequencies.
`nsave` | integer | -1 | Output cadence of restart dumps. If negative, no restart dumps are written.
`save_for_restart` | boolean | `false` | Write restart dumps. If `true`, then `nsave` should also be positive.
`write_phi_vs_time` | boolean | `false` | Write the full electrostatic potential \\( \varphi \\) to the output netCDF file.
`write_gvmus` | boolean | `false` | Writes \\(\sum_\boldsymbol{k}L_z^{-1}\int \textrm{d}z \lvert g_\boldsymbol{k} \rvert^2 \\) to the netCDF file. Resulting quantity varies over \\(v_\parallel\\), \\( \mu_s \\) and the species index \\( s\\).
`write_gzvs` | boolean | `false` |  Writes \\(\sum_\boldsymbol{k}L_z^{-1}\int \textrm{d} \mu_s \lvert g_\boldsymbol{k} \rvert^2 \\) to the netCDF file. Resulting quantity varies over \\(v_\parallel\\), \\( z \\) and the species index \\( s\\).
`write_omega` | boolean | `false` | Write the real and imaginary mode frequencies to both ASCII formats (`RUN_NAME.omega`) and the netCDF file.
`write_kspectra` | boolean | `false` | Writes \\(\langle \lvert \varphi(k_x, k_y, z)\rvert^2 \rangle_\psi \\) to the netCDF file.
`write_moments` | boolean | `false` | Write the moments of the distribution function (density, parallel velocity and temperature for each species) to the netCDF file.
`write_radial_fluxes` | boolean | `false` | Write the flux-surface-averaged radial fluxes to the netCDF file.
`write_radial_moments` | boolean | `false` | Write the flux-surface-averaged fluctuation moments to the netCDF file.
`write_fluxes_kxkyz` | boolean | `false` | Write the mode-by-mode radial fluxes as a function of \\( z \\) to the netCDF file.
`flux_norm` | boolean | `true` | If true, then scale radial fluxes by \\( \langle\lvert \nabla r\rvert \rangle_\psi \\), otherwise perform no rescaling. 
`nc_mult` | integer | 1 | Multiplication factor of the output cadence for netCDF files, which tend to take up more storage space than the ASCII output files.


# namelist `multibox_parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`boundary_size` | integer | 4 | Number of collocation points in the boundary region.
`krook_size` | integer | 0 | Number of collocation points in the Krook region. Will automatically max out at `boundary_size`.
`zf_option` | string | `'default'` | *Experimental* - Set how the zonal mode is handled during the communication of the multiple-flux-tube boundary condition. Should be one of <ul><li>`default` Do not treat the \\( k_y = 0\\) mode specially. </li><li>`skip_ky0` Do not communicate *any* information for the  \\( k_y = 0\\) mode, i.e., let it evolve without any input from auxilliary simulations. </li><li>`zero_ky0`Zero the \\( k_y = 0\\) mode in the boundary region.</li><li>`zero_fsa` Zero only the flux-surface-averaged component of the \\(k_y = 0\\) mode (keeping \\( v_{\parallel} \\) and \\( \mu_s \\) constant) in the boundary region. Note that this is ***not*** equivalent to zeroing the transit/bounce-averaged component, which would be the physically correct approach.</li></ul>
`krook_option` | string | `'exp'` | Shape of the Krook operator in the Krook region. Should be one of <ul><li>`flat` constant shape </li><li>`linear` linearly decreasing towards physical region </li><li>`exp` decreasing exponentially fast to the physical region, *recommended*</li><li>`exp_rev` decreasing exponentially slowly to the physical region.</li><li>`default` same as `exp`</li></ul>
`RK_step` | boolean | `false` | Communicate the multiple-flux-tube boundary condition at every implicit and Runge-Kutta substep, rather than one per timestep. 
`nu_krook_mb` | real | 0.0 | Strength (i.e. damping rate) of the Krook operator in the boundary region.
`mb_debug_step` | integer | -1 | If positive, print out gnuplot-readable binary dump of the real-space electrostatic potential \\( \varphi \)) at the outboard midplane at every `mb_debug_step` steps.
`krook_exponent` | real | 0.0 | Add a \\((k_y/k_{y0})^{\texttt{krook_exponent}}\\) prefactor to the multibox Krook operator.
`comm_at_init` | boolean | `false` | Communicate the boundary condition at the beginning of the simulation before time_advance is initialized.
`phi_bound` | integer | 0 | How many points into the boundary region on which to solve for the electrostatic potential \\( \varphi \\) when `ky_solve_real` is true.
`phi_pow` | real | 0.0 | *Experimental* - Which derivative of the electrostatic potential \\( \varphi \\) to solve for when `ky_solve_real` is true.
`krook_efold` | real | 3.0 | How many e-folds to use when `krook_option` is `exp` or `exp_rev`.
`use_dirichlet_BC` | boolean | `false` | Apply the Dirichlet boundary condition when performing a *single* flux-tube simulation, i.e., utilize the machinery in the multibox module, but use zeros in the communication instead of data generated from additional local simulations.
`LR_debug_option` | string | `'default'`| Left/Right debug option. Should be one of <ul><li>`default` Do not use this debug option.</li><li>`L` set \\( x \\) uniformly to \\( x_\mathrm{L} \\), the leftmost point in the physical region.</li><li>`R` set \\( x \\) uniformly to \\( x_\mathrm{R} \\), the rightmost point in the physical region.</li></ul> Useful for comparing the center domain to left/right domains in radially global simulation.
`smooth_ZFs` | boolean | `false` | *Experimental* - Try smoothing the zonal flows across the interface between boundary and physical region. *Never really worked properly.* 


# namelist `sources`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`source_option` | string | `'none'` | Type of source used for radially global simulation. Should be one of <ul><li> `none` no source used. </li><li>`default` same as `none`.</li><li>`Krook` Krook operator based source.</li><li>`projection` Projection operator based source.</li></ul> 
`conserve_momentum` | boolean | `false` | Enforce momentum conservation in the source.
`conserve_density` | boolean | `false` |  Enforce density conservation in the source.
`tcorr_source` | real  | 0.02 | Time correlation of the source. *Should be a longer timescale than any of interest in the problem.*
`nu_krook` | real | 0.05 | Strength (i.e. damping rate) of the Krook source. 
`ikxmax_source` | integer | 2 for `periodic_variation`, 1 otherwise | Maximum radial wavenumber on which the source acts.
`krook_odd` | boolean | `true` | *Experimental* - Only act on the modes that can affect the profiles.
`exclude_boundary_regions` | boolean | `radial_variation` and not `periodic_variation` | Ensure that source only acts on the physical region of the radial domain.
`tcorr_source_qn` | real | 0.0 | *Experimental* - Time correlation of the source in quasineutrality.
`exclude_boundary_regions_qn` | boolean | `exclude_boundary_regions` | Exclude boundary regions when applying the source in quasineutrality.
`from_zero` | boolean | `true` | Time correlated source starts from zero, rather than value of the distribution function at \\( t = 0 \\). 


# namelist `dissipation`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`include_collisions` | boolean | `true` | Include particle collisions.
`collisions_implicit` | boolean | `true` | Evaluate the collision operator implicitly.
`collision_model` | string | `'dougherty'` | Which collision operator to use. Options are `'dougherty'` (simplified operator) or `'fokker-planck'` (physical operator).
`hyper_dissipation` | boolean | `false` | Include hyper-dissipation. *Strongly* recommended for nonlinear simulations.


# namelist `collisions_dougherty`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`momentum_conservation` | boolean | `true` | Enforce momentum conservation.
`energy_conservation` | boolean | `true` | Enforce energy conservation.
`vpa_operator` | boolean | `true` | Include \\( \partial_{v_\parallel} \\) components of collision operator.
`mu_operator` | boolean | `true` | Include \\( \partial_\mu \\) components of collision operator.


# namelist `collisions_fp`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`testpart` | boolean | `true` | test particle component of Fokker-Planck operator. **Must be true.**
`fieldpart` | boolean | `false` | enable the field particle component (FPO) of the Fokker-Planck operator.
`interspec` | boolean | `true` | inter-species collisions in the Fokker-Planck operator.
`intraspec` | boolean | `true` | intra-species collisions in the Fokker-Planck operator.
`lmax` | integer | 1 | maximum l in spherical harmonic expansion of the field particle operator
`jmax` | integer | 1 | maximum j in Hirshman-Sigmar expansion of the field particle operator
`iiknob` | real | 1.0 | control the ion-ion collision frequency in Fokker-Planck operator
`ieknob` | real | 1.0 | control the ion-electron collision frequency in Fokker-Planck operator.
`eiknob` | real | 1.0 | control the electron-ion collision frequency in Fokker-Planck operator.
`eeknob` | real | 1.0 | control the electron-electron collision frequency in Fokker-Planck operator.
`eiediffknob` | real | 1.0 | control the electron-ion energy diffusion in Fokker-Planck operator.
`eideflknob ` | real | 1.0 | control the electron-ion pitch angle scattering in Fokker-Planck operator. 
`deflknob ` | real | 1.0 | control pitch angle scattering in Fokker-Planck operator, must be 1 or 0.
`eimassr_approx ` | boolean  | `false` | use mass ratio approxfimation for test particle operator, *beta*.
`advfield_coll ` | boolean | `true` | disable electrostatic potential terms in the field particle operator, *beta*.
`spitzer_problem ` | boolean | `false`| Solve the Spitzer problem for tests of the collision operator
`density_conservation` | boolean | `false` | if `True` and `equally_spaced_mu_grid=True` and `conservative_wgts_vpa=True`, then test-particle operator conserves density to machine precision.
`density_conservation_field` | boolean | `false` |  if `True` and `jmax`, `lmax < 2`, then field-particle operator conserves density to machine precision.
`density_conservation_tp` | boolean | `false` | if True add term to field particle operator to ensure density conservation, also on non-uniform grids.
`exact_conservation` | boolean | `false` | if True and `fieldpart=True` and `lmax=jmax=1` then momentum and energy conserved to machine precision - *in beta*. Works only if `nux = 0`, need to correct the discretisation of nux terms in test-particle operator.
`exact_conservation_tp` | boolean | `false` |  if True and `lmax=jmax=1` then momentum and energy conserved to machine precision, by using the test particle operator  to compute field particle terms; this is slower than exact_conservation.
`vpa_operator ` | boolean | `true`|  Include \\( \partial_{v_\parallel} \\) components of collision operator.
`mu_operator` | boolean | `true` | Include \\( \partial_\mu \\) components of collision operator.
`cfac` | real | 1.0 | scale gyrodiffusive term in test particle component of Fokker-Planck operator.
`cfac2` | real | 1.0 |  scale gyrodiffusive terms in field particle component of Fokker-Planck operator - *in beta*.
`nuxfac` | real | 1.0 | scale nux (mixed derivative) terms in test particle component of Fokker-Planck operator.
`i1fac` | real | 1.0 | for Spitzer problem
`i2fac` | real | 0.0 | for Spitzer problem
`no_j1l1` | boolean | `true` | disable j1l1 term in the field particle component of Fokker-Planck operator
`no_j1l2` | boolean | `false` | disable j1l2 term
`no_j0l2` | boolean | `false` | disable j0l2 term
`nvel_local` | integer | 512 | Size of velocity grid used for debugging. Currently, this option has no effect.


# namelist `hyper`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`D_hyper` | real | 0.05 | Maximal hyperdissipation damping rate.
`use_physical_ksqr` | boolean  | `true` if global, `false` otherwise | If true, use actual \\( k^2_\perp = k_x^2 \lvert \nabla x \rvert^2 + 2 k_xk_y (\nabla x \cdot \nabla y) + k_y^2 \lvert\nabla y \rvert^2 \\). Otherwise, use \\( k_\perp^2 = k_y^2[1 + (\theta - \theta_0)^2]\\).
`scale_to_outboard` | boolean | `false` | If true, scale maximal damping rate to maximum \\( k_\perp^2 \\) at outboard midplane. Otherwise, scale maximal damping to maximum \\( k_\perp^2 \\) over the entire domain.


# namelist `neoclassical_input`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`include_neoclassical_terms` | boolean | `false`| Include neoclassical terms. Need to use sfincs for input.
`neo_option` | string | `'sfincs'`| Option for obtaining neoclassical distribution function and potential. Currently, only `'sfincs'` is supported.
`nradii` | integer  | 5 | Number of radial points used for radial derivatives of neoclassical quantities.
`drho` | real | 0.01 | spacing in `rhoc` between points used for radial derivatives.


# namelist `sfincs_input`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`read_sfincs_output_from_file` | boolean | `false` | if true, will try to read in `Phi1Hat` and `delta_f` from pre-saved file named `sfincs.output`. Otherwise, run sfincs to compute these quantities on the fly.
`nproc_sfincs` | integer | 1 | Number of processors to use for sfincs calculations.
`calculate_radial_electric_field` | boolean | `true` | If true, then scan in radial electric field to find value for which ambipolarity is satisfied, and then use this value to obtain neoclassical fluxes, distribution function, and potential.
`includeXDotTerm` | boolean | `true` | include radial electric field term
`includeElectricFieldTermInXiDot` | boolean | `true` | 
`irad_min` | integer | `-nardii / 2` | minimum radial index (`irad=0` corresponds to central radius).
`irad_max` | integer | `nradii / 2` | maximum radial index.
`magneticDriftScheme` | integer | 0 | ???
`includePhi1` | boolean | `true` | If `true`, then `Phi1` will be calculated using quasineutrality.
`includePhi1InKineticEquation` | boolean | `false` |
`geometryScheme` | integer  | 1 | will be overridden by direct input of geometric quantities unless `geometryScheme = 5` (VMEC equilibrium)
`VMECRadialOption` | integer | 0 | **only relevant if `geometryScheme = 5`.** Radial option to use for VMEC equilibrium. Must be one of 
`equilibriumFile` | string | `'wout_161s1.nc'`| path of VMEC equilibrium file.
`coordinateSystem` | integer | 3 | ???
`inputRadialCoordinate` | integer | 3 | option 3 corresponds to using square root of toroidal flux normalized by toroidal flux enclose by the LCFS.
`inputRadialCoordinateForGradients` | integer | 3 | option 3 corresponds to same choise when calculating gradients of density, temperature, and potential
`ahat` | real | 1.0 | corresponds to \\( r_\mathrm{LCFS} \\) as reference length in sfincs. Only used in sfincs when `geometryScheme = 1`.
`psiAHat` | real | `geo_surf%psitor_lcfs` | \\( \psi_\mathrm{LCFS} / B_\mathrm{ref}a^2_\mathrm{ref} \\).
`Delta` | real | -1.0| \\( \Delta = \rho_\ast m_\mathrm{ref}v_\mathrm{th,ref}/(eB_\mathrm{ref}a_\mathrm{ref}) \\), with reference quantities given in SI units unless `geometryScheme` = 5, in which case \\(B_\mathrm{ref} = 1 \\) Tesla and \\(a_\mathrm{ref}  = 1 \\) meter (these are hardwired in sfincs). Set negative to allow check later to see if any value given in input file.
`dPhiHatdrN` | real | -9999.9| radial derivative of normalized \\( \phi \\).
`nu_N` | real | -1.0 | \\(\nu_n = \nu_\mathrm{ref} a_\mathrm{ref}/ v_\mathrm{th,ref}\\). Note that \\( \nu_\mathrm{ref} = 4 \sqrt{2 \pi} n_\mathrm{ref}e^4 \ln \Lambda / (3m^{1/2}_\mathrm{ref}T_\mathrm{ref}^{3/2})\\). Set negative to allow check later to see if any value given in input file.
`nxi` | integer | 48 | number of spectral coefficients in pitch angle
`nx` | integer  | 12 | number of speeds
`ntheta` | integer | 65 | number of poloidal angles. 
`nzeta` | integer | 1 | number of toroidal angles. 1 is approriate for tokamaks.
`Er_window` | real | 0.3 | ???


