!###############################################################################
!#################### READ STELLA NAMELISTS FOR NEOCLASSICS ####################
!###############################################################################
! 
! This module will read the namelists associated with neoclassical terms:
! 
!    neoclassical_input
!       include_neoclassical_terms = .false.
!       neo_option = 'sfincs'
!       nradii = 5.0
!       drho = 0.01
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
module namelist_neoclassical_input

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_neoclassical_input

   ! Parameters need to be public
   public :: neo_option_sfincs

   private

   ! Create parameters for <neo_option>
   integer, parameter :: neo_option_sfincs = 1

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist
   
contains

   !****************************************************************************
   !                             NEOCLASSICAL INPUT                            !
   !****************************************************************************
   subroutine read_namelist_neoclassical_input(include_neoclassical_terms, neo_option_switch, nradii, drho)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: neo_option_switch
      integer, intent (out) :: nradii
      real, intent (out) :: drho
      logical, intent(out) :: include_neoclassical_terms

      ! Local variable to set <neo_option_sfincs>
      character(10) :: neo_option
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_neoclassical_input
      call read_input_file_neoclassical_input
      call check_inputs_neoclassical_input

   contains

      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_neoclassical_input

         implicit none

         include_neoclassical_terms = .false.   ! set to .true. to include neoclassical terms in GK equation
         nradii = 5.0                           ! number of radial points used for radial derivative ! of neoclassical quantities
         drho = 0.01                            ! spacing in rhoc between points used for radial derivatives
         neo_option = 'sfincs'                  ! option for obtaining neoclassical distribution function and potential

      end subroutine set_default_parameters_neoclassical_input

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_neoclassical_input

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr
         
         ! Link text options for <neo_option> to an integer value
         type(text_option), dimension(2), parameter :: neoopts = (/ &
             text_option('default', neo_option_sfincs), &
             text_option('sfincs', neo_option_sfincs)/)

         ! Variables in the <neoclassical_input> namelist
         namelist /neoclassical_input/ include_neoclassical_terms, neo_option, nradii, drho
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('neoclassical_input', dexist)
         if (dexist) read (unit=in_file, nml=neoclassical_input)

         ! Read the text option in <neo_option> and store it in <neo_option_switch>
         ierr = error_unit()
         call get_option_value(neo_option, neoopts, neo_option_switch, &
             ierr, 'neo_option in neoclassical_input')

      end subroutine read_input_file_neoclassical_input

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_neoclassical_input

         implicit none

         if (nradii /= 3 .and. nradii /= 5) then
             write (*, *) 'WARNING: only nradii of 3 or 5 is currently supported in neoclassical_input namelist'
             write (*, *) 'WARNING: forcing nradii=5'
             nradii = 5
         end if

      end subroutine check_inputs_neoclassical_input

   end subroutine read_namelist_neoclassical_input

end module namelist_neoclassical_input
