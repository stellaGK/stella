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
