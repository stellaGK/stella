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
   <bmag> = B / B_ref
   <gradx_dot_gradx> = |∇x|²
   <grady_dot_grady> = |∇y|²
   <gradx_dot_grady> = ∇x . ∇y
   <B_times_gradB_dot_gradx> = B × ∇B · ∇x (a*B_ref/B^3)
   <B_times_gradB_dot_grady> = B × ∇B · ∇y (a*B_ref/B^3)
   <B_times_kappa_dot_gradx> = B × κ · ∇x (a*B_ref/B^2)
   <B_times_kappa_dot_grady> = B × κ · ∇y (a*B_ref/B^2)
   <gradpar> = a nabla_parallel z

The normalized derivatives are defined as,
   <dxdpsi> = (rho_r/a) (d tilde{x} / d tilde{psi})
   <dydalpha> = (rho_r/a) (d tilde{y} / d tilde{alpha})
