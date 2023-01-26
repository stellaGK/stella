# namelist `dissipation`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`include_collisions` | boolean | `true` | Include particle collisions.
`collisions_implicit` | boolean | `true` | Evaluate the collision operator implicitly.
`collision_model` | string | `'dougherty'` | Which collision operator to use. Options are `'dougherty'` (simplified operator) or `'fokker-planck'` (physical operator).
`hyper_dissipation` | boolean | `false` | Include hyper-dissipation. *Strongly* recommended for nonlinear simulations.
