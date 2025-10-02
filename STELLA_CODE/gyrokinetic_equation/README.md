# Gyrokinetic equation

Stella evolves the normalized gyrokinetic equation in time, give by equation
(19) in [2019 - Barnes - stella an operator-split, implicit-explicit df-gyrokinetic 
code for general magnetic field configurations], where each variable is assumed to be normalized:

$$ 
\begin{aligned}
\frac{\partial \tilde{g}_{\mathbf{k},s}}{\partial \tilde{t}} 
{+} \tilde{v}_{\text{th},s} \tilde{v}_{\parallel}\, \hat{b} \cdot \tilde{\boldsymbol{\nabla}} \tilde{z} \, \Big( \frac{\partial \tilde{g}_{\mathbf{k},s}}{\partial \tilde{z}} + \frac{Z_s}{\tilde{T}_s} \frac{\partial J_0(a_{\mathbf{k},s}) \tilde{\varphi}_\mathbf{k}}{\partial \tilde{z}} \, e^{-\tilde{v}^2_s} \Big)
{-} \tilde{v}_{\text{th},s}  \tilde{\mu} \, \hat{b} \cdot \tilde{\boldsymbol{\nabla}} \tilde{B} \,\frac{\partial \tilde{g}_{\mathbf{k},s}}{\partial \tilde{v}_\parallel} \\ 
{+} i \omega_{d,\mathbf{k},s} \Big(  \tilde{g}_{\mathbf{k},s} 
{+} \frac{Z_s}{\tilde{T}_s} J_0(a_{\mathbf{k},s}) \tilde{\varphi}_\mathbf{k} e^{-\tilde{v}^2_s} \Big) 
{+} i \omega_{*,\mathbf{k},s} J_0(a_{\mathbf{k},s}) \tilde{\varphi}_\mathbf{k} 
{+} \mathcal{N}_{\mathbf{k},s} = 0
\end{aligned}
$$

## Terms in the gyrokinetic equation

In the gyrokinetic equation given above, there are 6 terms
   - (1) Time derivative
   - (2) Parallel streaming because of b . ∇z
   - (3) Magnetic mirror because of b . ∇B
   - (4) Magnetic drift because of omega_{d,k,s}
   - (5) Driving term, since omega_{*,k,s} depends on the density and temperature gradients
   - (6) The nonlinear term since different (kx,ky) modes can interact through this term
   
### Parallel streaming

The parallel streaming term is evolved in 'gk_parallel_streaming.f90'.
   
   
### Magnetic mirror

The magnetic mirror term is evolved in 'gk_mirror.f90'.

As ions and electrons move freely along magnetic field lines, their magnetic moment 
(mu = v_perp²/2B) and kinetic energy (m v_perp²/2 + m v_parallel²/2) are conserved.
As a result, if the magnetic field strength B(z) increases along z, the perpendicular
velocity increases in order to preserve the magnetic moment, and the parallel velocity
of the particle decreases to conserve the kinetic energy. As a result, it is possible
that there is a point along the magnetic field line, where B(z) is such that the parallel
velocity of the particle becomes zero, and the particle reverses direction. Therefore,
these so-called trapped particles bounce back and forward between two bounce points, 
also called mirror points. Hence, the term containing b . ∇B is the mirror term.

   
### Magnetic drift

The drift term is evolved in 'gk_drifts.f90'.


### Nonlinear term

The nonlinear term is evolved in 'gk_nonlinearity.f90'.

It is common to turn off the nonlinear term to study the behaviour of individual
(kx,ky) modes in so-called linear simulations.


## Time advance gyrokinetic equation

The gyrokinetic equation is evolved in time in 'gk_time_advance.f90'.


## Write documentation for:
   
   - response_matrix.fpp
   - timestep_calculations.f90
   - gk_radial_variation.f90
   - gk_sources.fpp
   - gk_ffs_solve.f90
   - gk_flow_shear.f90
   - gk_implicit_solve.f90
   
   
## Clebsch form

We assume the equations are derived using
    B = C ∇ψ x ∇α

If the radial coordinate is chosen to be psi = psi_t
    B = sign_torflux ∇ψ x ∇α 
    <clebsch_factor> = sign_torflux
    <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha)
    <drhodpsi> = drho/dψ̃ = d(r/a)/d(psi/(a^2*Br)) = (a*Bref) * dr/dpsi

If the radial coordinate is chosen to be psi = psi_p
    B = - ∇ψ x ∇α
    <clebsch_factor> = -1
    <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha) 
    <drhodpsi> = drho/dψ̃ = d(r/a)/d(psi/(a^2*Br)) = (a*Bref) * dr/dpsi

If the radial coordinate is chosen to be psi = q 
    B = - (dpsi_p/dq) ∇ψ x ∇α 
    1/<clebsch_factor> = - dq/d(psip/(a^2*Br)) = - (a^2*Bref) (qd/dpsi_p)
    <drhodpsi> = drho/dq = d(r/a)/dq = (1/a) * dr/dq
    <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha) 
