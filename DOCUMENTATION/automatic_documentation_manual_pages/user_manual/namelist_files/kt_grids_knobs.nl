# namelist `kt_grids_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`grid_option` | string | `'default'`| Sets the layout of the perpendicular grid. Should be one of <ul><li>`range` use a set range of wavenumbers. Suitable for linear simulation.</li> <li>`box` use a physical box in coordinate space, i.e., a wavenumber grid that satisfies \\(A_\boldsymbol{k} = A^\ast_{-\boldsymbol{k}} \\). *Required* for nonlinear simulation.</li><li>`default` same as `range`</li><li>`annulus` same as `box`</li><li>`nonlinear` same as `box`</li></ul>
