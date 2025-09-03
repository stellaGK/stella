!###############################################################################
!################### READ STELLA NAMELISTS FOR VELOCITY GRID ###################
!###############################################################################
! 
! This module will read the namelists associated with the velocity grid:
! 
!   velocity_grids
!     nvgrid = 24.0
!    nmu = 12.0
!    vpa_max = 3.0
!    vperp_max = 3.0
!    equally_spaced_mu_grid = .false.
!    conservative_wgts_vpa = .false.
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
!###############################################################################
module namelist_velocity_grids

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_velocity_grids

   private

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                               VELOCITY GRIDS                              !
   !****************************************************************************
   subroutine read_namelist_velocity_grids (nvgrid, nmu, vpa_max, vperp_max, &
      equally_spaced_mu_grid, conservative_wgts_vpa)

      use mp, only: proc0

      implicit none 

      ! Variables that are read from the input file
      integer, intent(out) :: nvgrid
      integer, intent(out) :: nmu
      real, intent(out) :: vpa_max, vperp_max
      logical, intent(out) :: equally_spaced_mu_grid, conservative_wgts_vpa
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_velocity_grids
      call read_input_file_velocity_grids
      call write_parameters_to_input_file

   contains
      !------------------------ Default input parameters -----------------------
      subroutine set_default_velocity_grids

         implicit none

         nvgrid = 24
         vpa_max = 3.0
         nmu = 12
         vperp_max = 3.0
         equally_spaced_mu_grid = .false.
         conservative_wgts_vpa = .false.

      end subroutine set_default_velocity_grids

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_velocity_grids

         use file_utils, only: input_unit_exist

         implicit none

         ! Variables in the <velocity_grids> namelist
         namelist /velocity_grids/ nvgrid, nmu, vpa_max, vperp_max, &
               equally_spaced_mu_grid, conservative_wgts_vpa
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('velocity_grids', dexist)
         if (dexist) read (unit=in_file, nml=velocity_grids)

      end subroutine read_input_file_velocity_grids
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&velocity_grids'
         write (unit, '(A, I0)') '  nvgrid = ', nvgrid
         write (unit, '(A, I0)') '  nmu = ', nmu
         write (unit, '(A, F0.2)') '  vpa_max = ', vpa_max
         write (unit, '(A, F0.2)') '  vperp_max = ', vperp_max
         write (unit, '(A, L0)') '  equally_spaced_mu_grid = ', equally_spaced_mu_grid
         write (unit, '(A, L0)') '  conservative_wgts_vpa = ', conservative_wgts_vpa
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_velocity_grids

end module namelist_velocity_grids
