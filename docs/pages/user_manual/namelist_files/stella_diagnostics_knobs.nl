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
