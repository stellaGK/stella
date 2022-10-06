# namelist `parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`beta`  | float | 0.0 | Plasma \\( \beta \\). Currently has no effect.
`zeff`  | float | 1.0 | Effective charge number for use with *effective* electron-ion and electron-impurity collisions in the Fokker-Planck collision operator (see `ecoll_zeff`).
`tite`  | float | 1.0 | Ratio of ion to electron temperature, \\( T_\mathrm{i}/T_\mathrm{e} \\). Used in quasineutrality when adiabatic species is used.
`nine`  | float |1.0 | Ratio of ion to electron density, \\( n_\mathrm{i}/n_\mathrm{e} \\). Used in quasineutrality when adiabatic species is used.
`rhostar`  | real | -1.0 | The gyrokinetic expansion parameter \\( \rho_\mathrm{th,ref}/a_\mathrm{ref} \\). For effects beyond the flux-tube limit (full-flux-surface, radially global, neoclassical terms, etc...). Overwritten if `irhostar` is positive.
`irhostar`  | real | -1.0 | Sets `rhostar = 1.0 / irhostar` if positive.
`vnew_ref`  | real | -1.0 | Reference collision frequency. Various input options will overwrite this if it is negative.
`g_exb`  | real | 1.0 | Equilibrium \\( \boldsymbol{E \times B} \\) shear rate. More specifically, \\( \gamma_\boldsymbol{ E \times B} = (r/q) (\textrm{d}\omega / \textrm{d}r) R_0/\sqrt{2}v_\mathrm{th,ref}\\). Uses the Hammett wavenumber shift method, with nonlinear corrections proposed by McMillan.
`g_exbfac`  | real | 1.0 | Prefactor for perpendicular component of equilibrium \\( \boldsymbol{E \times B} \\) flow shear. Setting to 0.0 turns this component off.
`omprimfac`  | real | 1.0 | Prefactor for parallel component of equilibrium \\( \boldsymbol{E \times B} \\) flow shear. Setting to 0.0 turns this component off.
