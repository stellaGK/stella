# namelist `vpamu_grid_parameters`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`nvgrid` | integer | 24 | Number of positive points in \\( v_\parallel \\). Actual grid will be twice as large. Note that \\( v_\parallel = 0\\) is not included.
`nmu` | integer | 12 | Number of points in the magnetic moment \\( \mu \\).
`vpa_max` | real | 3.0 | Maximum \\( v_\parallel \\) on the grid in terms of thermal velocity.
`vperp_max` | real | 3.0 | Maximum \\( v_\perp \\) on the grid in terms of thermal velocity.
`equally_spaced_mu_grid` | boolean | `false` | If false, use Gaussian quadrature points for \\( \mu \\) grid (recommended). Otherwise, use equally spaced grid. Both options exclude \\( \mu = 0 \\).
`conservative_wgts_vpa` | boolean | `false` | Use density-conserving weights for \\( v_\parallel \\). Used primarily by Fokker-Planck collision operator.
