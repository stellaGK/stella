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
