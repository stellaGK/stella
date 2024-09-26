# namelist `species_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`nspec` | integer | 2  | Number of species.
`species_option` | string | `'stella'` | How to read in species data. Should be one of <ul><li>`stella` Read from stella input file.</li><li>  `euterpe` Read form euterpe file. </li><li> `input.profiles` Reads in General Atomics `input.gacode` file. **This may need to be updated for newer `gacode` files.**</li><li>`default` same as `stella`</li></ul> 
`read_profile_variation` | boolean | false | Save information necessary for recalculating kinetic profile quantities away from \\( r \\). Used for performing local simulations at different radial locations in order to compare to global simulation.
`write_profile_variation` | boolean | false | Recomputes the kinetic profile information using kinetic profile quantities stored in a file, originally calculated at some \\( r_\mathrm{file} \\), at new \\( r \\).
`ecoll_zeff` | boolean | false | If true, use an effective intra-species electron-ion collision rate using `zeff`.
