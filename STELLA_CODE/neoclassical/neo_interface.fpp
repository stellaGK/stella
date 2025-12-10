! ################################################################################# !
! ################################ NEO INTERFACE ################################## !
! ################################################################################# !
!
! A module for reading in NEO output data, to be used as input data for second order
! gyrokinetic runs in stella. 
!
! NEO output on the distribution correction, H_1, and the electrostatic potential
! correction, ϕ^1_0, is related to F_1 via:
!
!
!       F_1 = H_1 - ( e * Z / T ) * ϕ^1_0 * F_0
!
! where:
!   - F_1: Total first-order correction to the distribution function.
!   - H_1: NEO's distribution function correction.
!   - ϕ^1_0: NEO's electrostatic potential correction.
!   - F_0: The Maxwellian background distribution.
!   - T, Z, e: Temperature, charge number, and elementary charge.
!
! stella requires three seperate NEO runs (for three seperate flux surfaces) 
! to calculate the gradient drive term arising from F_1. Files to be read in are: 
!  - out.neo.f
!  - out.neo.f.right
!  - out.neo.f.left
!  - out.neo.phi
!  - out.neo.phi.right
!  - out.neo.phi.left
!  - out.neo.grid 
!  - out.neo.equil
!  - out.neo.species
!  - out.neo.version
!
! ################################################################################# !

module neo_interface



end module neo_interface 
