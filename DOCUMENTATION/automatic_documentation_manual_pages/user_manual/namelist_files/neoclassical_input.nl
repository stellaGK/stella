# namelist `neoclassical_input`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`include_neoclassical_terms` | boolean | `false`| Include neoclassical terms. Need to use sfincs for input.
`neo_option` | string | `'sfincs'`| Option for obtaining neoclassical distribution function and potential. Currently, only `'sfincs'` is supported.
`nradii` | integer  | 5 | Number of radial points used for radial derivatives of neoclassical quantities.
`drho` | real | 0.01 | spacing in `rhoc` between points used for radial derivatives.
