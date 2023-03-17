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
