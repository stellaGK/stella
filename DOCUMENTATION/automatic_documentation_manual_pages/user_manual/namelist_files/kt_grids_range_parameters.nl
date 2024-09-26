# namelist `kt_grids_range_parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`naky` | integer | 1  | Number of points along the binormal mode grid direction.
`nakx` | integer | 1  | Number of points along the radial mode grid direction.
`aky_min` | real | 0.0  | Minimum binormal wavenumber, or binormal wavenumber if `naky` is 1.
`aky_max` | real | 0.0 | Maximum binormal wavenumber if `naky` larger than 1.
`akx_min` | real | 0.0  | Minimum radial wavenumber. Unused if \\( \theta_0 \\) is properly specified.
`akx_max` | real | -1.0 | Maximum radial waveumner.
`theta0_min` | real | 0.0  | Minimum \\( \theta_0 = k_x / \hat{s} k_y \\). Unused if \\( \theta_{0,\textrm{ min}} > \theta_{0,\textrm{ max}} \\).
`theta0_max` | real | -1.0  | Maximum \\( \theta_0 \\).
`kyspacing_option` | string | `'default'`  | Sets the spacing between \\(a k_y\\) grid points, available options are : <ul><li>`default` same as `linear`</li><li>`linear` linear spacing in \\(a k_y\\)</li><li>`exponential` linear spacing in \\(\log{(a k_y)}\\)</ul>
