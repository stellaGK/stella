!###############################################################################
!                                    ARRAYS                                    
!###############################################################################
! This module stores any useful arrays which are used throughout the code.
!###############################################################################
module arrays

   use mpi
   use stella_common_types, only: response_matrix_type

   implicit none
   
   ! Velocity-dependent coefficients used in the equations
   public :: kperp2, dkperp2dr
   public :: time_gke
   public :: time_parallel_nl
   public :: initialised_radial_variation
   public :: initialised_wdrift, initialised_wstar
   public :: initialised_parallel_streaming, initialised_implicit_drifts
   public :: wstar, wstarp
   public :: wdriftx_g, wdrifty_g
   public :: wdriftx_phi, wdrifty_phi
   public :: wdriftx_bpar, wdrifty_bpar
   public :: wdriftpx_g, wdriftpy_g
   public :: wdriftpx_phi, wdriftpy_phi

   ! Arrays without velocity dependence. Used mostly in field calculations.
   public :: response_matrix, response_window
   public :: shift_state
   public :: denominator_QN, ddenominator_QNdr
   public :: denominator_QN13, denominator_QN_MBR1, denominator_QN_MBR3
   public :: denominator_QNinv11, denominator_QNinv13, denominator_QNinv31, denominator_QNinv33
   public :: denominator_QN_MBR
   public :: apar_denom
   public :: theta
   public :: c_mat
   public :: exclude_boundary_regions_qn
   public :: tcorr_source_qn, exp_fac_qn
   public :: qn_window, qn_zf_window
   public :: denominator_QN_h, denominator_QN_MBR_h, efac, efacp
   public :: time_field_solve
   
   !----------------------------------------------------------------------------
   ! For the gyrokinetic equation
   !----------------------------------------------------------------------------
   
   ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)
   real, dimension(:, :, :), allocatable :: wstar, wstarp

   ! (nalpha, -nzgrid:nzgrid, -vmu-layout-)
   real, dimension(:, :, :), allocatable :: wdriftx_g, wdrifty_g
   real, dimension(:, :, :), allocatable :: wdriftx_phi, wdrifty_phi
   real, dimension(:, :, :), allocatable :: wdriftx_bpar, wdrifty_bpar 
   real, dimension(:, :, :), allocatable :: wdriftpx_g, wdriftpy_g
   real, dimension(:, :, :), allocatable :: wdriftpx_phi, wdriftpy_phi

   ! dkperp2dr will contain the radial variation of kperp2
   real, dimension(:, :, :, :), allocatable :: kperp2, dkperp2dr

   ! for time advance
   real, dimension(2, 10) :: time_gke = 0.
   real, dimension(2, 2) :: time_parallel_nl = 0.

   logical :: initialised_wdrift
   logical :: initialised_wstar
   logical :: initialised_parallel_streaming
   logical :: initialised_radial_variation
   logical :: initialised_implicit_drifts

   !----------------------------------------------------------------------------
   ! For field solves
   !----------------------------------------------------------------------------
   
   type(response_matrix_type), dimension(:), allocatable :: response_matrix
   integer :: response_window = MPI_WIN_NULL
   real, dimension(:), allocatable :: shift_state
   
   ! The electrostatic potential phi is calculated based on the quasi-neutrality condition
   !     sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g + (Zs/Ts) (Gamma0 - 1) phi ] = 0
   !     phi = sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J0 * g ] / [ sum_s (Zs²ns/Ts) (1 - Gamma0) ]
   !     denominator_QN[iky,ikz,iz] = sum_s (Zs²ns/Ts) (1 - Gamma0)
   ! The denominators needed to calculate <phi> are initialised in 
   !     - quasineutrality_equation_fluxtube::init_quasineutrality_equation_fluxtube
   real, dimension(:, :, :), allocatable :: denominator_QN, ddenominator_QNdr
   real, dimension(:, :, :), allocatable :: denominator_QN13, denominator_QN_MBR1, denominator_QN_MBR3
   real, dimension(:, :, :), allocatable :: denominator_QNinv11, denominator_QNinv13, denominator_QNinv31, denominator_QNinv33
   real, dimension(:, :), allocatable :: denominator_QN_MBR
   real, dimension(:, :, :), allocatable ::  apar_denom
   
   ! (nakx, nakx, -nzgrid:nzgrid)
   complex, dimension(:, :, :), allocatable :: theta
   
   ! (nakx, nakx)
   complex, dimension(:, :), allocatable :: c_mat
   
   ! Variables needed for the sources
   logical :: exclude_boundary_regions_qn
   real :: tcorr_source_qn, exp_fac_qn
   integer :: qn_window = MPI_WIN_NULL, qn_zf_window = MPI_WIN_NULL
   real :: denominator_QN_h, denominator_QN_MBR_h, efac, efacp
   real, dimension(2, 5) :: time_field_solve = 0.

end module arrays
