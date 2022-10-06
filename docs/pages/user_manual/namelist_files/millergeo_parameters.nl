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
