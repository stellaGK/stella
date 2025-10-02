!###############################################################################
!                                    TIMERS                                    
!###############################################################################
! This module stores the timers for the routines.
!###############################################################################
module timers

   implicit none
   
   public :: time_total
   public :: time_init
   public :: time_gke
   public :: time_mirror
   public :: time_sources
   public :: time_multibox
   public :: time_collisions
   public :: time_parallel_nl
   public :: time_field_solve
   public :: time_all_diagnostics
   public :: time_implicit_advance
   public :: time_parallel_streaming
   public :: time_individual_diagnostics
   public :: time_response_matrix
   public :: time_lu_decomposition
   
   private
   
   !----------------------------------------------------------------------------

   ! Time the entire stella code
   real, dimension(2) :: time_total = 0.
   
   ! Time the initialisation of stella
   real, dimension(2) :: time_init = 0.
   
   ! Time the response matrix
   real, dimension(2) :: time_response_matrix = 0.
   real, dimension(2) :: time_lu_decomposition = 0.
   
   ! Time the gyrokinetic equation
   !    - time_gke(:,1) = Gyrokinetic equation + field equations
   !    - time_gke(:,2) = unused
   !    - time_gke(:,3) = unused
   !    - time_gke(:,4) = Magnetic drift (omega_d) wdrifty
   !    - time_gke(:,5) = Magnetic drift (omega_d) wdriftx
   !    - time_gke(:,6) = Drive term (omega_*)
   !    - time_gke(:,7) = ExB nonlinear term
   !    - time_gke(:,8) = Explicit terms
   !    - time_gke(:,9) = Implicit terms
   !    - time_gke(:,10) = Radial variation
   real, dimension(2, 10) :: time_gke = 0.
   
   ! The the mirror term
   !    - time_mirror(:,1) = Mirror term
   !    - time_mirror(:,2) = Redistribute
   real, dimension(2, 2) :: time_mirror = 0.
   
   ! The the parallel nonlinearity
   !    - time_parallel_nl(:,1) = Parallel nonlinearity
   !    - time_parallel_nl(:,2) = Redistribute
   real, dimension(2, 2) :: time_parallel_nl = 0.
   
   ! The the collisions
   !    - time_collisions(:,1) = Collisions
   !    - time_collisions(:,2) = Redistribute
   real, dimension(2, 2) :: time_collisions = 0.
   
   ! Time the implicit advance
   !    - time_implicit_advance(:,1) = Implicit time advance
   !    - time_implicit_advance(:,2) = Bidiagonal solve
   !    - time_implicit_advance(:,3) = Back substitution
   real, dimension(2, 3) :: time_implicit_advance = 0.
   real, dimension(2, 3) :: time_parallel_streaming = 0.

   ! Time the field equations routines
   !    - time_field_solve(:,1) = Evolve the fields
   !    - time_field_solve(:,2) = Redistribute
   !    - time_field_solve(:,3) = int_dv_g
   !    - time_field_solve(:,4) = calculate_phi
   !    - time_field_solve(:,5) = calculate_phi_adia_elec
   real, dimension(2, 5) :: time_field_solve = 0.
   
   ! Time the diagnostics
   !    - time_all_diagnostics(:) = all diagnostics
   !    - time_individual_diagnostics(:,1) = omega
   !    - time_individual_diagnostics(:,2) = potential
   !    - time_individual_diagnostics(:,3) = omega
   !    - time_individual_diagnostics(:,4) = fluxes
   !    - time_individual_diagnostics(:,5) = moments
   !    - time_individual_diagnostics(:,6) = distribution
   real, dimension(2) :: time_all_diagnostics = 0.
   real, dimension(2, 6) :: time_individual_diagnostics = 0.
   
   ! Time sources
   !    - time_sources(:,1) = sources
   !    - time_sources(:,2) = redistribute
   real, dimension(2, 2) :: time_sources = 0.
   
   ! Time multi box communications
   !    - time_multibox(:,1) = mb_comm
   !    - time_multibox(:,2) = mb_krook
   real, dimension(2, 2) :: time_multibox = 0.

end module timers
