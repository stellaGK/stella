# namelist `layouts_knobs`

These options control the ordering of the coordinate and velocity variables when the domain is decomposed. **At the moment, there is no reason to use anything other than the default options.**

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`xyzs_layout`  | string | `"yxzs"` | Order of the coordinate variables when decomposed for multicore simulations.
`vms_layout` | string | `"vms"` | Order of the velocity variables when decomposed for multicore simulations.
