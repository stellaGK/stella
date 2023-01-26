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
