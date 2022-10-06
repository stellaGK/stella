# namelist `sources`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`source_option` | string | `'none'` | Type of source used for radially global simulation. Should be one of <ul><li> `none` no source used. </li><li>`default` same as `none`.</li><li>`Krook` Krook operator based source.</li><li>`projection` Projection operator based source.</li></ul> 
`conserve_momentum` | boolean | `false` | Enforce momentum conservation in the source.
`conserve_density` | boolean | `false` |  Enforce density conservation in the source.
`tcorr_source` | real  | 0.02 | Time correlation of the source. *Should be a longer timescale than any of interest in the problem.*
`nu_krook` | real | 0.05 | Strength (i.e. damping rate) of the Krook source. 
`ikxmax_source` | integer | 2 for `periodic_variation`, 1 otherwise | Maximum radial wavenumber on which the source acts.
`krook_odd` | boolean | `true` | *Experimental* - Only act on the modes that can affect the profiles.
`exclude_boundary_regions` | boolean | `radial_variation` and not `periodic_variation` | Ensure that source only acts on the physical region of the radial domain.
`tcorr_source_qn` | real | 0.0 | *Experimental* - Time correlation of the source in quasineutrality.
`exclude_boundary_regions_qn` | boolean | `exclude_boundary_regions` | Exclude boundary regions when applying the source in quasineutrality.
`from_zero` | boolean | `true` | Time correlated source starts from zero, rather than value of the distribution function at \\( t = 0 \\). 
