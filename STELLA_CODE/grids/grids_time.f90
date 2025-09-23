!###############################################################################
!                              TIME GRID MODULE                                                 
!###############################################################################
! This module initilaises dt for stella. 
!###############################################################################m
module grids_time

   implicit none

   public :: code_dt, update_time, code_dt_old
   public :: code_time
   public :: write_dt
   public :: init_tstart, init_delt, checkcodedt
   public :: cfl_dt_linear, cfl_dt_ExB, cfl_dt_parallel
   public :: code_dt_min, code_dt_max

   private

   real :: cfl_dt_linear = -1.
   real :: cfl_dt_ExB = -1.
   real :: cfl_dt_parallel = -1.
   real :: code_dt, code_dt_min, code_dt_max

   ! added May 18, 2009 to take care of problems
   ! in exb_shear calculation after change in time step size
   real :: code_dt_old = 0.
   real :: code_time = 0.

contains

   !****************************************************************************
   !                         Initialise start time
   !****************************************************************************
   subroutine init_tstart(tstart)

      real, intent(in) :: tstart

      code_time = tstart

   end subroutine init_tstart

   !****************************************************************************
   !                        Initialise time step size
   !****************************************************************************
   subroutine init_delt(delt, delt_max, delt_min)
      real, intent(in) :: delt, delt_max, delt_min

      code_dt = delt
      code_dt_min = delt_min

      ! Do not allow code_dt to increase beyond the input value
      ! Unless we specified delt_max in the input file
      ! For example, when restarting the simulation
      if (delt_max < 0) then
         code_dt_max = code_dt
      else
         code_dt_max = delt_max
      end if

   end subroutine init_delt

   subroutine update_time
! MAB+CMR, 21/5/09: set code_dt_old to code_dt BEFORE any changes in timestep
      code_dt_old = code_dt
      code_time = code_time + code_dt

   end subroutine update_time

   !****************************************************************************
   !                        Write time step information
   !****************************************************************************
   subroutine write_dt

      if (cfl_dt_linear > 0. .and. cfl_dt_linear < 1.e7) &
         write (*, *) 'TIME STEP:'
      write (*, '(A12, ES10.2E2)') "   cfl_dt:"//repeat(' ', 50), cfl_dt_linear
      write (*, '(A12, ES10.2E2)') "   code_dt:"//repeat(' ', 50), code_dt

   end subroutine write_dt

   !****************************************************************************
   !                     Check that code_dt is not too small
   !****************************************************************************
   subroutine checkcodedt(stop_stella)

      use mp, only: proc0, broadcast
      logical, intent(in out) :: stop_stella

      if (proc0) then
         if (code_dt < code_dt_min) then
            stop_stella = .true.
         end if
      end if

      if (proc0 .and. (code_dt < code_dt_min)) then
         write (*, *)
         write (*, *) 'EXITING STELLA BECAUSE CODE_DT<CODE_DT_MIN:'
         write (*, '(A16, ES10.2E2)') "   code_dt:"//REPEAT(' ', 50), code_dt
         write (*, '(A16, ES10.2E2)') "   code_dt_min:"//REPEAT(' ', 50), code_dt_min
      end if

      call broadcast(stop_stella)

   end subroutine checkcodedt

end module grids_time
