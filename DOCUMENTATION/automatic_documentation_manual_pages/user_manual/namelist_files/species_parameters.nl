# namelist `species_parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`z` | integer  | 1 | charge number.
`mass` | real | 1.0 | particle mass, normalized to \\( m_\mathrm{ref} \\). 
`dens` | real | 1.0 | species density, normalized to \\( n_\mathrm{ref} \\). 
`temp` | real | 1.0 | species temperature, normalized to \\( T_\mathrm{ref} \\). 
`fprim` | real | -999.9 | Density gradient scale length \\(a_\mathrm{ref}L_{n_s}^{-1} = - a_\mathrm{ref} \textrm{d} \ln n_s/ \textrm{d}r \\)). Note the negative sign.
`tprim` | real | -999.9 | Temperature gradient scale length \\(a_\mathrm{ref}L_{T_s}^{-1} = - a_\mathrm{ref} \textrm{d} \ln T_s/ \textrm{d}r \\)). Note the negative sign.
`d2ndr2` | real | 0.0 | Second derivative of density, normalized to \\( a_\mathrm{ref}^2/n_\mathrm{ref}\\). For use with radially global simulation.
`d2Tdr2` | real | 0.0 | Second derivative of temperature, normalized to \\( a_\mathrm{ref}^2/T_\mathrm{ref}\\). For use with radially global simulation.
`bess_fac` | real | 1.0 | Perfactor for Bessel function argument. Setting to 0.0 renders particle drift-kinetic.
`type` | string | `'default'` | Particle type. Should be one of <ul><li>`ion` ion species.</li><li>`default` same as `ion`. </li><li>  `electron` electron species. </li><li> `e` same as `electron` </li><li>`beam` slowing down species. </li><li>`fast` same as `beam`. </li><li>`alpha` same as `beam`. </li><li>`slowing-down` same as `beam`. </li><li>`trace` tracer species. </li></ul> 
