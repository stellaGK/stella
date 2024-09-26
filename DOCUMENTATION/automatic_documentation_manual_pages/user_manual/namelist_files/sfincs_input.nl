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
