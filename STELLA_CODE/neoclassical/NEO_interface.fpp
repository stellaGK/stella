! ##################################################################################################### !
! ########################################## NEO INTERFACE ############################################ !
! ##################################################################################################### !
!
! This module reads in NEO output data associated with the first order neoclassical correction to the 
! equilibrium distributuion, F_1. This is needed as input for second order simulations in stella.   
! 
! NEO output for the distribution correction, H_1, and the electrostatic potential 
! correction, ϕ^1_0, is related to F_1 via
!
! F_1 = H_1 - (e * Z/T) * ϕ^1_0 * F_0
!
! where:
!   - F_1: Total first-order correction to the distribution function.
!   - H_1: NEO's distribution function correction.
!   - ϕ^1_0: NEO's electrostatic potential correction.
!   - F_0: The Maxwellian background distribution.
!   - e, Z, T: Elementary Charge, charge number, and temperature.
!
! stella requires output data from 3 seperate NEO runs (for 3 neighbouring flux surfaces). This is needed
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
! Files without a "left" or "right" suffix correspond to the central flux surface for which stella
! simulations are to be evaluated.
!
! ##################################################################################################### !

module NEO_interface

    implicit none

    public :: read_basic_neo_files, read_neo_f_and_phi, neo_version_data, neo_equil_data
    public :: neo_grid_data, neo_species_data

    private


! ##################################################################################################### !
! ################# Represents the contents of out.neo.version - NEO metadata. ######################## ! 
! ##################################################################################################### !

    type neo_version_data
        character(len=:), allocatable :: commit
        character(len=:), allocatable :: system
        character(len=:), allocatable :: date
    end type neo_version_data


! ##################################################################################################### !
! ################## Represents the contents of out.neo.grid - NEO grid data. ######################### ! 
! ##################################################################################################### !

    type neo_grid_data
        integer :: n_species = -1
        integer :: n_energy = -1
        integer :: n_xi = -1
        integer :: n_theta = -1
        real(8), dimension(:), allocatable :: theta
        integer :: n_radial = -1
        real(8), dimension(:), allocatable :: radius
  end type neo_grid_data


! ##################################################################################################### !
! ################ Represents the contents of out.neo.equil - NEO equilibrium data. ################### !
! ##################################################################################################### !

  type neo_equil_data
     real(8), dimension(:), allocatable :: radius  
     real(8), dimension(:), allocatable :: radial_electric_field
     real(8), dimension(:), allocatable :: q_safety
     real(8), dimension(:), allocatable :: rho_star
     real(8), dimension(:), allocatable :: major_radius
     real(8), dimension(:), allocatable :: angular_frequency
     real(8), dimension(:), allocatable :: rotation_shear ! 1-dimensional arrays for radially varying quantities. 
     
     real(8), dimension(:, :), allocatable :: density
     real(8), dimension(:, :), allocatable :: temperature
     real(8), dimension(:, :), allocatable :: density_gradient
     real(8), dimension(:, :), allocatable :: temperature_gradient
     real(8), dimension(:, :), allocatable :: collision_frequency ! 2-D arrays for radius and species dependent quantites. 
  end type neo_equil_data


! ##################################################################################################### !
! ################# Represents the contents of out.neo.species - NEO species data. #################### !
! ##################################################################################################### !

  type neo_species_data
     real(8), dimension(:), allocatable :: mass
     real(8), dimension(:), allocatable :: charge
  end type neo_species_data


contains

end module NEO_interface



