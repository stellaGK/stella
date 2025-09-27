!###############################################################################
!                                    TIMERS                                    
!###############################################################################
! This module stores the timers for the routines.
!###############################################################################
module timers

   implicit none
   
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
   
   private
   
   !----------------------------------------------------------------------------

   ! Time the gyrokinetic equation routines
   real, dimension(2, 10) :: time_gke = 0.
   real, dimension(2, 2) :: time_mirror = 0.
   real, dimension(2, 2) :: time_collisions = 0.
   real, dimension(2, 2) :: time_parallel_nl = 0.
   real, dimension(2, 3) :: time_implicit_advance = 0.
   real, dimension(2, 3) :: time_parallel_streaming = 0.

   ! Time the field equations routines
   real, dimension(2, 5) :: time_field_solve = 0.
   
   ! Time the diagnostics 
   real, dimension(2) :: time_all_diagnostics = 0.
   real, dimension(2, 6) :: time_individual_diagnostics = 0.
   
   ! Time radial variation and sources
   real, dimension(2, 2) :: time_sources = 0.
   real, dimension(2, 2) :: time_multibox = 0.

end module timers
