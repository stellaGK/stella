# Geometry

This module construct the geometric quantities needed in stella. When investigating
the gyrokinetic equation, it is clear that we need to define b · ∇z and B as a
function of (psi, alpha, z). Moreover, within the drift frequency the terms
B × ∇B · ∇x, B × ∇B · ∇y, B × κ · ∇x and B × κ · ∇y appear. Finally in the norm
of the perpendicular wavenumber, which is needed for the Bessel functions, the terms
|∇x|², |∇y|² and ∇x . ∇y appear, all of which are a function of (psi, alpha, z).

In the flux-tube approximation, a specific field line on a specific flux surface
are chosen, hence we have a fixed (psi0, alpha0) and the geometric quantities
are only a function of z, which is the coordinate parallel to the magnetic field 
lines. Common choices for z include:
    - The poloidal angle, θ, for Miller equilibria
    - The toroidal angle, ζ, for VMEC equilibria
    - The arc length along the field line, ℓ


## Geometric quantities

The following geometric quantities are defined in this module:
   <bmag>(alpha,z) = B / B_ref
   <gradx_dot_gradx>(alpha,z) = |∇x|²
   <grady_dot_grady>(alpha,z) = |∇y|²
   <gradx_dot_grady>(alpha,z) = ∇x . ∇y
   <B_times_gradB_dot_gradx>(alpha,z) = B × ∇B · ∇x (a*B_ref/B^3)
   <B_times_gradB_dot_grady>(alpha,z) = B × ∇B · ∇y (a*B_ref/B^3)
   <B_times_kappa_dot_gradx>(alpha,z) = B × κ · ∇x (a*B_ref/B^2)
   <B_times_kappa_dot_grady>(alpha,z) = B × κ · ∇y (a*B_ref/B^2)
   <b_dot_gradz>(alpha,z) = b · ∇z
   <b_dot_gradz_avg>(z) = sum_alpha b · ∇z J dalpha / sum_alpha J dalpha

The normalized derivatives are defined as,
   <dxdpsi> = (rho_r/a) (d tilde{x} / d tilde{psi})
   <dydalpha> = (rho_r/a) (d tilde{y} / d tilde{alpha})
