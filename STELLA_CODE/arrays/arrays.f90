!###############################################################################
!                                    ARRAYS                                    
!###############################################################################
! This module stores any useful arrays which are used throughout the code.
! These are often global arrays, or arrays that are used in multiple modules,
! so they are stored here to make them easily accessible.
!###############################################################################
module arrays
   
   ! Import variables for the response matrix
   use common_types, only: response_matrix_type
   use mpi, only: MPI_WIN_NULL

   implicit none
   
   ! Keep track of which routines have been initialised
   public :: initialised_wdrift
   public :: initialised_wstar
   public :: initialised_wpol
   public :: initialised_parallel_streaming
   public :: initialised_implicit_drifts
   public :: initialised_radial_variation
   
   !----------------------------------------------------------------------------
   ! For the Gyrokinetic Equation
   !----------------------------------------------------------------------------
   
   ! Velocity-dependent ferquencies used in the equations
   public :: wstar
   public :: wstarp
   public :: wpol
   public :: wdriftx_g
   public :: wdrifty_g
   public :: wdriftx_phi
   public :: wdrifty_phi
   public :: wdriftx_bpar
   public :: wdrifty_bpar
   public :: wdriftpx_g
   public :: wdriftpy_g
   public :: wdriftpx_phi
   public :: wdriftpy_phi
   
   ! Geometric quantities vs (naky, nakx, nalpha, -nzgrid:nzgrid)
   public :: kperp2
   public :: dkperp2dr

   ! For radial variations and sources
   public :: theta
   public :: c_mat
   public :: exclude_boundary_regions_qn
   public :: tcorr_source_qn
   public :: exp_fac_qn
   
   ! For flow shear
   public :: shift_state
   
   ! For the reponse matrix
   public :: response_matrix
   public :: response_window
   public :: qn_window
   public :: qn_zf_window
   
   !----------------------------------------------------------------------------
   ! For the Field Equations
   !----------------------------------------------------------------------------
   
   ! Arrays used to calculate the fields for electrostatic simulations
   public :: denominator_fields
   public :: denominator_fields_h
   
   ! Arrays used to calculate the fields for electrostatic simulations
   ! considering a Modified Boltzmann Response for the electrons
   public :: denominator_fields_MBR
   public :: denominator_fields_MBR_h
   public :: efac, efacp
   
   ! Arrays used to calculate the fields for electromagnetic simulations
   public :: denominator_fields_inv11
   public :: denominator_fields_inv13
   public :: denominator_fields_inv31
   public :: denominator_fields_inv33
   public :: apar_denom
   
   ! Arrays used to calculate the fields for radial variation simulations
   public :: denominator_fields_dr
   
   private
   
   !----------------------------------------------------------------------------

   ! Keep track of which routines have been initialised
   logical :: initialised_wdrift
   logical :: initialised_wstar
   logical :: initialised_wpol
   logical :: initialised_parallel_streaming
   logical :: initialised_radial_variation
   logical :: initialised_implicit_drifts
   
   !----------------------------------------------------------------------------
   ! For the Gyrokinetic Equation
   !----------------------------------------------------------------------------
   
   ! Frequencies appearing in the gyrokinetic equations vs (nalpha, -nzgrid:nzgrid, -vmu-layout-)
   real, dimension(:, :, :), allocatable :: wstar, wstarp
   real, dimension(:, :, :), allocatable :: wpol
   real, dimension(:, :, :), allocatable :: wdriftx_g, wdrifty_g
   real, dimension(:, :, :), allocatable :: wdriftx_phi, wdrifty_phi
   real, dimension(:, :, :), allocatable :: wdriftx_bpar, wdrifty_bpar
   real, dimension(:, :, :), allocatable :: wdriftpx_g, wdriftpy_g
   real, dimension(:, :, :), allocatable :: wdriftpx_phi, wdriftpy_phi

   ! Geometric quantities vs (naky, nakx, nalpha, -nzgrid:nzgrid)
   real, dimension(:, :, :, :), allocatable :: kperp2
   real, dimension(:, :, :, :), allocatable :: dkperp2dr
   
   ! For the reponse matrix
   type(response_matrix_type), dimension(:), allocatable :: response_matrix
   integer :: qn_window = MPI_WIN_NULL
   integer :: qn_zf_window = MPI_WIN_NULL
   integer :: response_window = MPI_WIN_NULL
   
   ! Variables needed for the sources and radial variation
   complex, dimension(:, :, :), allocatable :: theta ! (nakx, nakx, -nzgrid:nzgrid)
   complex, dimension(:, :), allocatable :: c_mat    ! (nakx, nakx)
   logical :: exclude_boundary_regions_qn
   real :: tcorr_source_qn, exp_fac_qn
   
   ! For flow shear
   real, dimension(:), allocatable :: shift_state

   !----------------------------------------------------------------------------
   ! For the Field Equations
   !----------------------------------------------------------------------------
   
   ! For electrostatic simulations
   real, dimension(:, :, :), allocatable :: denominator_fields     ! (nakx, naky, -nzgrid:nzgrid)
   real, dimension(:, :), allocatable :: denominator_fields_MBR    ! (nakx, -nzgrid:nzgrid)
   real :: denominator_fields_h
   real :: denominator_fields_MBR_h
   real :: efac, efacp
   
   ! For electromagnetic simulations (nakx, naky, -nzgrid:nzgrid)
   real, dimension(:, :, :), allocatable :: denominator_fields_inv11
   real, dimension(:, :, :), allocatable :: denominator_fields_inv13
   real, dimension(:, :, :), allocatable :: denominator_fields_inv31
   real, dimension(:, :, :), allocatable :: denominator_fields_inv33
   real, dimension(:, :, :), allocatable :: apar_denom
   
   ! For radial variation simulations (nakx, naky, -nzgrid:nzgrid)
   real, dimension(:, :, :), allocatable :: denominator_fields_dr

end module arrays
