# Gyrokinetic equation

Stella evolves the normalized gyrokinetic equation in time, give by equation
(19) in [2019 - Barnes - stella an operator-split, implicit-explicit df-gyrokinetic 
code for general magnetic field configurations], where each variable is assumed to be normalized:

   (1)      ∂g_{k,s} / ∂t =
   (2)      - v_{th,s} v_parallel b . ∇z ( ∂g_{k,s} / ∂z + Z_s/T_s ∂(J_0 ϕ_k) / ∂z exp(-v*2_s) 
   (3)      + v_{th,s} mu_s b . ∇B ∂g_{k,s} / ∂v_parallel
   (4)      - i omega_{d,k,s} (g_{k,s} + Z_s/T_s J_0 ϕ_k exp(-v*2_s) 
   (5)      - i omega_{*,k,s} J_0 ϕ_k 
   (6)      - (B_r/2) (dy/dalpha) (dx/dpsi) F_k ( F^{-1}_k [ik_y J_0 ϕ_k] F^{-1}_k [ik_x g_{k,s}] - F^{-1}_k [ik_x J_0 ϕ_k] F^{-1}_k [ik_y g_{k,s} ] )
            
If maxwellian_normalization = .true., the evolved distribution function is normalised by a Maxwellian.

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

! We assume the equations are derived using
!     B = C ∇ψ x ∇α
! 
! If the radial coordinate is chosen to be psi = psi_t
!     B = sign_torflux ∇ψ x ∇α 
!     <clebsch_factor> = sign_torflux
!     <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha)
!     <drhodpsi> = drho/dψ̃ = d(r/a)/d(psi/(a^2*Br)) = (a*Bref) * dr/dpsi
! 
! If the radial coordinate is chosen to be psi = psi_p
!     B = - ∇ψ x ∇α
!     <clebsch_factor> = -1
!     <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha) 
!     <drhodpsi> = drho/dψ̃ = d(r/a)/d(psi/(a^2*Br)) = (a*Bref) * dr/dpsi
! 
! If the radial coordinate is chosen to be psi = q 
!     B = - (dpsi_p/dq) ∇ψ x ∇α 
!     1/<clebsch_factor> = - dq/d(psip/(a^2*Br)) = - (a^2*Bref) (qd/dpsi_p)
!     <drhodpsi> = drho/dq = d(r/a)/dq = (1/a) * dr/dq
!     <dydalpha> = (rhor/a)(d(y/rhor)/dalpha) = (1/a)(dy/dalpha) 
