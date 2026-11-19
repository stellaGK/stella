!###############################################################################
!############ READ STELLA NAMELISTS FOR RADIAL VARIATION (MULTIBOX) ############
!###############################################################################
! 
! This module will read the namelists associated with shear flow:
!   
!   flow_shear
!     prp_shear_enabled = .false.
!     hammett_flow_shear = .true.
!     g_exb = 0.0
!     g_exbfac = 1.0
!     omprimfac = 1.0
!     remap_em = 0
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
! 
!###############################################################################
module namelist_flow_shear

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_flow_shear

   private

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   
   !****************************************************************************
   !                               FLOW SHEAR TERMS                            !
   !****************************************************************************
   subroutine read_namelist_flow_shear(prp_shear_enabled, hammett_flow_shear, &
      g_exb, g_exbfac, omprimfac, remap_em)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      logical, intent(out) :: prp_shear_enabled, hammett_flow_shear 
      integer, intent(out) :: remap_em
      real :: g_exb, g_exbfac, omprimfac
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_flow_shear
      call read_input_file_flow_shear
      call check_inputs_flow_shear
      call write_parameters_to_input_file

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_flow_shear

         implicit none

         prp_shear_enabled = .false.
         hammett_flow_shear = .true.
         g_exb = 0.0
         g_exbfac = 1.0
         omprimfac = 1.0
         remap_em = 0

      end subroutine set_default_parameters_flow_shear

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_flow_shear

         use mp, only: mp_abort

         implicit none

         if (remap_em /= 0 .and. remap_em /= 1) then
            call mp_abort('flow_shear: remap_em must be 0 (remap g) or 1 (remap gbar). Aborting.')
         end if

      end subroutine check_inputs_flow_shear

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_flow_shear

         use file_utils, only: input_unit_exist
         implicit none

         namelist /flow_shear/ prp_shear_enabled, hammett_flow_shear, g_exb, g_exbfac, omprimfac, remap_em
         in_file = input_unit_exist('flow_shear', dexist)
         if (dexist) read (unit=in_file, nml=flow_shear)

      end subroutine read_input_file_flow_shear
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&flow_shear'
         write (unit, '(A, L1)') '  prp_shear_enabled = ', prp_shear_enabled
         write (unit, '(A, L1)') '  hammett_flow_shear = ', hammett_flow_shear
         write (unit, '(A, ES0.4)') '  g_exb = ', g_exb
         write (unit, '(A, ES0.4)') '  g_exbfac = ', g_exbfac
         write (unit, '(A, ES0.4)') '  omprimfac = ', omprimfac
         write (unit, '(A, I0)') '  remap_em = ', remap_em
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_flow_shear
   
end module namelist_flow_shear
