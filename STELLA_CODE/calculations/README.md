# Calculations

This directory contains common calculations performed in stella.
This includes, derivative calculations (both finite difference, and spectral),
gyro-averages, and volume averages.

## calculations_add_explicit_terms

This module adds explicit terms in the time advance routine. It is used frequently,
therefore, this calculation is written once and called upon for (e.g.)

  - advance_wdriftx_explicit
  - advance_wdrifty_explicit
  - advance_wstar_explicit

## calculations_checksum

This module calculates the sum of the modulus squared of the field and
distribution function arrays. This is helpful for debugging and ensuring
that the arrays are being correctly updated.

## calculations_finite_differences

This module is used for computing derivatives using finite difference schemes.
There are different schemes available, including upwind schemes of various orders.

## calculations_gyro_averages

This module computes the gyro-average of the fields and distribution funcitons.

## calculations_kxky_derivatives

This module computes spectral derivatives in kx, ky.

## calculations_kxky

This module swaps between different orderings of the (kx, ky) grids.

Use reality to swap between arrays with:
    ky >= 0 and all kx (ordered like 0, ..., kxmax, -kxmax, ..., -dkx)
    kx >= 0 and all ky (ordered like 0, ..., kymax, -kymax, ..., -dky)

This is needed because Fourier transforms may need one or the other grid.
If transforming from ky->y then we need all ky, and if transforming from
kx->x then we need all kx.

## calculations_redistribute

This module swaps between different distributions:

    kxkyz2vmu - swap between (ky, kx, z) local with (vpa, mu, species) parallelised
                and (vpa, mu) local with (ky, kx, z, species) parallelised. This is also
                used for the reverse swap (i.e. velocity local to spatial coordinate local)
     kxyz2vmu - swap between (y, kx, z) local with (vpa, mu, species) parallelised
                and (vpa, mu) local with (y, kx, z, species) parallelised. This is also
                used for the reverse swap (i.e. velocity local to spatial coordinate local).
                Note, this is needed for full flux annulus as y often needs to be held in real space.
      xyz2vmu - swap between (x, y, z) local with (vpa, mu, species) parallelised and
                (vpa, mu) local with (x, y, z) parallised. This is also used for the reverse swap.
   kymus2vmus - swap between (kx, z, vpa) local with (ky, mu, species) parallelised and (ky, kx, z)
                local with (vpa, mu, species) parallelised. This is also used for the reverse swap.

## calculations_timestep

This module computes CFL conditions on the time step size, and resets it when needed.

## calculations_tofrom_ghf

This module swaps between different forms of the distribution funciton.

    f - The distribution function evaluated at particle position.
    g - The gyroaveraged distribution funciton, evaluated at gyro-center position.
    h - The non-adiabatic part of the distribution function, evaluated at gyro-center position.

Take r to be particle position, and R to be the gyro-center position. Then, for the electrostatic case:

    g(R) = <f(r)>_R, where <.>_R defines the gyro-average at fixed gyro-center
    f(r) = h(R) - Z_s/T_s * phi(r) * F_0 
    g(R) = h(R) - Z_s/T_s * <phi(r)>_R * F_0 = f(r) + Z_s/T_s * [phi(r) - <phi(r)>_R] * F_0 

We need these, because although stella mostly works with g(R), some calclations require h or f.

## calculations_transforms

This module performs discrete Fourier transformations of both real and complex data.

## calculations_velocity_integrals

This module computes any integration in velocity space,

    int d^3 v f = int_0^{2*pi} d \sigma int_{-infty}^{infty} d vpa \int_0^{infty} d mu B/m f

## calculations_volume_averages

This module compute spatial integrations, such as field line averages and volume integrals.
