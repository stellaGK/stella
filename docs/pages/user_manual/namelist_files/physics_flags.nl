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
