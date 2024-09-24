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
