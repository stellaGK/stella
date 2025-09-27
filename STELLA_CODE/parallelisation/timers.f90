!###############################################################################
!                                    TIMERS                                    
!###############################################################################
! This module stores the timers for the routines.
!###############################################################################
module timers

   implicit none
   
   public :: time_gke
   public :: time_parallel_nl
   public :: time_field_solve
   
   private
   
   !----------------------------------------------------------------------------

   ! Time the gyrokinetic equation routines
   real, dimension(2, 10) :: time_gke = 0.
   real, dimension(2, 2) :: time_parallel_nl = 0.

   ! Time the field equations routines
   real, dimension(2, 5) :: time_field_solve = 0.

end module timers
