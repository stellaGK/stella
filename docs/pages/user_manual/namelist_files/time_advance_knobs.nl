# namelist `time_advance_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`xdriftknob` | real | 1.0 | Prefactor for radial magnetic drift. Setting to 0.0 turns the term off.
`ydriftknob` | real | 1.0 | Prefactor for the binormal magnetic drift. Setting to 0.0 turns the term off.
`wstarknob` | real | 1.0 | Prefactor for the \\( \omega_\ast \\) term. Setting to 0.0 turns the term off.
`explicit_option` | string | `'default'` | Chooses the Runge-Kutta scheme for the explicit integration. Should be one of <ul><li>`rk2` second-order Runge-Kutta. </li> <li>`rk3` third-order strong-stability-preserving Runge-Kutta (recommended). </li><li>`rk4` fourth-order Runge-Kutta. </li><li>`default` same as `rk3`. </li></ul> Note that higher-order Runge-Kutta schemes can increase memory usage.
`flip_flop` | boolean | `false` | Utilize the flip-flopping approach that flips the integration order every time-step. Should increase time accuracy, at least linearly. Does sometimes lead to spurious oscillations.
