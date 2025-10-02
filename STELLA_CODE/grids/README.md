# Grids

This document describes the grids used in stella, which define the basis variables of the system.
Stella uses the following grid structure:
            (ky, kx, z, vpa, mu, species) 

## Spatial Grids 
Stella works in {x,y,z} coordinate, which form a left-handed, orthonormal basis.

### ky, kx Grids
Here, (ky, kx) are the Fourier component of the y and x coordinates, respectively: 
            y = dy/dα (α - α_0)
            x = dx/dψ (ψ - ψ_0)
where (α_0, ψ_0) are the values of (α, ψ) at center of the simulation domain. 
Here, (α, ψ) are the field line choosing, and radial coordinates respectively. 
Both (ky, kx) are normalised by 1/ρ_r, where ρ_r is the reference gyroradius. 

There is some freedom in choosing x and y. The exact definitions of x and y depend 
on the chosen magnetic geometry. For details, refer to the geometry directory.

Regardless of the choice of geometry, (ky, kx) are chosen to be either a uniformly spaced,
or exponentially spaced grid in Fourier space. How these are initialised depends on if considering we are 'box' or 'range' mode.  Typically, 'range' is used for linear simulations, and 'box' is 
used for nonlinear simulations. 

Range mode takes in minimum and maximum values of (ky, kx) and initalises a uniform
grid in this range of values, depending on the number of modes specified for each, 
(N_ym N_x).

Box mode takes in values of y_0 and x_0. These are defined as:
            y_0 = 2 * pi * N_y * Δy
            x_0 = 2 * pi * N_x * Δx

Hence, with (y_0, x_0) and (N_y, N_x) specified, (Δy, Δx) can be computed. 
Then, (Δk_y, Δk_x) can be computed from:
            Δk_y = 1 / y_0
            Δk_x = 1 / x_0


### z grid

z is the field-aligned coordinate, measuring the distance along the magnetic field. The parallel
dynamics are treated using the real space coordinate z to correctly capture the parallel
boundary conditions.
Common choices for z are the toroidal angle, ζ, the poloidal angle, θ, and the arc length along 
the field line.


## Velocity Grids
For the velocity grids, stella uses (vpa, μ, σ). Here, σ is the gyration angle, which we average over
in the gyrokinetic formalism, so is not explicitly present in the computation. Hence, stella has grids only for
(vpa, μ).

### vpa grid

vpa is the velocity, projected along the magentic field:
            vpa = v · b̂

vpa is normalised by v_th in stella.

### μ grid

Here, μ is the magnetic moment, defined as: 
            μ = m v_{\perp}^2 / 2B

μ is normalised by 2*T_s / B_r, where T_s is the species temperature, and B_r is the reference magentic field. 


## Species grid

The species grid indexes the particle species (e.g., ions, electrons, or additional species).
Species can be defined in the stella input file or imported from an EUTERPE file.