! ##################################################################################################### !
! ########################################## NEO INTERFACE ############################################ !
! ##################################################################################################### !
!
! This module reads in NEO output data associated with the first order neoclassical correction to the 
! equilibrium distributuion, F_1. This is needed as input for second order simulations in stella.   
! 
! NEO output on the distribution correction, H_1, and the electrostatic potential 
! correction, ϕ^1_0, is related to F_1 via
!
! F_1 = H_1 - (e * Z/T) * ϕ^1_0 * F_0
!
! where:
!   - F_1: Total first-order correction to the distribution function.
!   - H_1: NEO's distribution function correction.
!   - ϕ^1_0: NEO's electrostatic potential correction.
!   - F_0: The Maxwellian background distribution.
!   - T, Z, e: Temperature, charge number, and elementary charge.
!
! stella requires output data from 3 seperate NEO runs for 3 neighbouring flux surfaces. This is needed
! to calculate the equilibrium gradient drive arising from F_1. Files to be read in are: 
!
! out.neo.f
! out.neo.f.right
! out.neo.f.left
! out.neo.phi
! out.neo.phi.right
! out.neo.phi.left
! out.neo.version
! out.neo.species
! out.neo.grid
! out.neo.version
!
! ##################################################################################################### !


