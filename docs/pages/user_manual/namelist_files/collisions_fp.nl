# namelist `collisions_fp`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`testpart` | boolean | `true` | test particle component of Fokker-Planck operator. **Must be true.**
`fieldpart` | boolean | `false` | enable the field particle component (FPO) of the Fokker-Planck operator.
`interspec` | boolean | `true` | inter-species collisions in the Fokker-Planck operator.
`intraspec` | boolean | `true` | intra-species collisions in the Fokker-Planck operator.
`lmax` | integer | 1 | maximum l in spherical harmonic expansion of the field particle operator
`jmax` | integer | 1 | maximum j in Hirshman-Sigmar expansion of the field particle operator
`iiknob` | real | 1.0 | control the ion-ion collision frequency in Fokker-Planck operator
`ieknob` | real | 1.0 | control the ion-electron collision frequency in Fokker-Planck operator.
`eiknob` | real | 1.0 | control the electron-ion collision frequency in Fokker-Planck operator.
`eeknob` | real | 1.0 | control the electron-electron collision frequency in Fokker-Planck operator.
`eiediffknob` | real | 1.0 | control the electron-ion energy diffusion in Fokker-Planck operator.
`eideflknob ` | real | 1.0 | control the electron-ion pitch angle scattering in Fokker-Planck operator. 
`deflknob ` | real | 1.0 | control pitch angle scattering in Fokker-Planck operator, must be 1 or 0.
`eimassr_approx ` | boolean  | `false` | use mass ratio approxfimation for test particle operator, *beta*.
`advfield_coll ` | boolean | `true` | disable electrostatic potential terms in the field particle operator, *beta*.
`spitzer_problem ` | boolean | `false`| Solve the Spitzer problem for tests of the collision operator
`density_conservation` | boolean | `false` | if `True` and `equally_spaced_mu_grid=True` and `conservative_wgts_vpa=True`, then test-particle operator conserves density to machine precision.
`density_conservation_field` | boolean | `false` |  if `True` and `jmax`, `lmax < 2`, then field-particle operator conserves density to machine precision.
`density_conservation_tp` | boolean | `false` | if True add term to field particle operator to ensure density conservation, also on non-uniform grids.
`exact_conservation` | boolean | `false` | if True and `fieldpart=True` and `lmax=jmax=1` then momentum and energy conserved to machine precision - *in beta*. Works only if `nux = 0`, need to correct the discretisation of nux terms in test-particle operator.
`exact_conservation_tp` | boolean | `false` |  if True and `lmax=jmax=1` then momentum and energy conserved to machine precision, by using the test particle operator  to compute field particle terms; this is slower than exact_conservation.
`vpa_operator ` | boolean | `true`|  Include \\( \partial_{v_\parallel} \\) components of collision operator.
`mu_operator` | boolean | `true` | Include \\( \partial_\mu \\) components of collision operator.
`cfac` | real | 1.0 | scale gyrodiffusive term in test particle component of Fokker-Planck operator.
`cfac2` | real | 1.0 |  scale gyrodiffusive terms in field particle component of Fokker-Planck operator - *in beta*.
`nuxfac` | real | 1.0 | scale nux (mixed derivative) terms in test particle component of Fokker-Planck operator.
`i1fac` | real | 1.0 | for Spitzer problem
`i2fac` | real | 0.0 | for Spitzer problem
`no_j1l1` | boolean | `true` | disable j1l1 term in the field particle component of Fokker-Planck operator
`no_j1l2` | boolean | `false` | disable j1l2 term
`no_j0l2` | boolean | `false` | disable j0l2 term
`nvel_local` | integer | 512 | Size of velocity grid used for debugging. Currently, this option has no effect.
