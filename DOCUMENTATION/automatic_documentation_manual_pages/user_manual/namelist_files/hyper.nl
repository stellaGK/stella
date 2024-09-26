# namelist `hyper`


Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`D_hyper` | real | 0.05 | Maximal hyperdissipation damping rate.
`use_physical_ksqr` | boolean  | `true` if global, `false` otherwise | If true, use actual \\( k^2_\perp = k_x^2 \lvert \nabla x \rvert^2 + 2 k_xk_y (\nabla x \cdot \nabla y) + k_y^2 \lvert\nabla y \rvert^2 \\). Otherwise, use \\( k_\perp^2 = k_y^2[1 + (\theta - \theta_0)^2]\\).
`scale_to_outboard` | boolean | `false` | If true, scale maximal damping rate to maximum \\( k_\perp^2 \\) at outboard midplane. Otherwise, scale maximal damping to maximum \\( k_\perp^2 \\) over the entire domain.
