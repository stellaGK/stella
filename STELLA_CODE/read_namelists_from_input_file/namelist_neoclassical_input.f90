! ======================================================================================================= !
! --------------------------- Read Stella Namelist for Neoclassical Physics. ---------------------------- !
! ======================================================================================================= !
! 
! This module will read the namelists associated with neoclassical terms:
! 
!    &neoclassical_input
!       include_neoclassical_terms = .false. ! < ==== Boolean for controlling the inclusion of  neoclassical corrections.
!       neo_option = 'NEO' ! String for controlling the neoclassical information provider: sfincs or NEO. 
!       nradii = 5 ! Integer number of flux surfaces in the drift kinetic calculations. Needed to calulcate F_1 gradients for the higher order equilibrium drive. Typically taken as !  or 5. 
!       drho = 0.01 ! Radial spacing between the flux surfaces. 
!    /
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
! First the default input parameters are set, then the default options are overwritten with those specified 
! in the input file. Optionally, it is checked whether any input variables are clashing with each other.
! 
! ======================================================================================================= !


module namelist_neoclassical_input
   implicit none

   ! Make reading routines accesible to other modules.
   public :: read_namelist_neoclassical_input

   ! Parameters need to be public.
   public :: neo_option_sfincs
   public :: neo_option_NEO

   private

   ! Create parameters for <neo_option>.
   integer, parameter :: neo_option_sfincs = 1
   integer, parameter :: neo_option_NEO = 2    ! <=========================== NEW PARAMETER FOR NEO OPTION!


   ! These variables are used in every single subroutine, so make them global.
   integer :: in_file
   logical :: dexist
   
contains


! ======================================================================================================= !
! ------------------------------------------ Neoclassical Input. ---------------------------------------- !
! ======================================================================================================= !

   subroutine read_namelist_neoclassical_input(include_neoclassical_terms, neo_option_switch, nradii, drho)

      use mp, only: proc0             ! Only using the master processor for reading. 

      implicit none

      ! Variables that are read from the input file. 
      integer, intent (out) :: neo_option_switch      ! Should either be 1 or 2 for sfincs or NEO, respectively.
      integer, intent (out) :: nradii
      real, intent (out) :: drho
      logical, intent(out) :: include_neoclassical_terms

      ! Local variable to set <neo_option_sfincs>.
      character(10) :: neo_option

      if (.not. proc0) return
      call set_default_parameters_neoclassical_input
      call read_input_file_neoclassical_input
      call check_inputs_neoclassical_input
      call write_parameters_to_input_file

   contains
   
! ======================================================================================================= !
! --------------------------------------- Default Input Parameters. ------------------------------------- !
! ======================================================================================================= !

      subroutine set_default_parameters_neoclassical_input
         implicit none

         include_neoclassical_terms = .false.   ! Set to .true. to include neoclassical terms in GK equation.
         nradii = 3                             ! Number of radial points used for radial derivative of neoclassical quantities. Typically taken as 3, may break for other integers. 
         drho = 0.01                            ! Spacing in rhoc between points used for radial derivatives.
         neo_option = 'NEO'                     ! Option for obtaining neoclassical distribution function and potential.

      end subroutine set_default_parameters_neoclassical_input


! ======================================================================================================= !
! --------------------------------------- Default Input Parameters. ------------------------------------- !
! ======================================================================================================= !

      subroutine read_input_file_neoclassical_input

         use file_utils, only: input_unit_exist
         use file_units, only: unit_error_file
         use text_options, only: text_option, get_option_value

         implicit none
         
         ! Link text options for <neo_option> to an integer value
         type(text_option), dimension(3), parameter :: neoopts = (/ &
             text_option('default', neo_option_NEO), &
             text_option('NEO', neo_option_NEO), &
             text_option('sfincs', neo_option_sfincs)/)

         ! Variables in the <neoclassical_input> namelist
         namelist /neoclassical_input/ include_neoclassical_terms, neo_option, nradii, drho

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('neoclassical_input', dexist)
         if (dexist) read (unit=in_file, nml=neoclassical_input)

         ! Read the text option in <neo_option> and store it in <neo_option_switch>
         call get_option_value(neo_option, neoopts, neo_option_switch, &
             unit_error_file, 'neo_option in neoclassical_input')    ! Takes the string neo_option from the namelist and matches it against the neoopts array. 
                                                                     ! This produces the integer neo_option_switch. Should be either 1 (for sfincs) or 2 (for NEO). 
                                                                     ! Other options are not supported. 

      end subroutine read_input_file_neoclassical_input


! ======================================================================================================= !
! ---------------------------------------- Check Input Parameters. -------------------------------------- !
! ======================================================================================================= !

      subroutine check_inputs_neoclassical_input
      
         use mp, only: mp_abort

         implicit none
         
         ! if (include_neoclassical_terms) then
         !    call mp_abort("Many mistakes are present in the definitions of gds23 and gds25. Aborting.")  <======== Hardcoded abort when sfincs option chosen? 
         ! end if 

         if (nradii /= 3) then
             write (*, *) 'WARNING: only nradii of 3 is currently supported in neoclassical_input namelist'
             write (*, *) 'WARNING: forcing nradii=3'
             nradii = 3
         end if

      end subroutine check_inputs_neoclassical_input


! ======================================================================================================= !
! ---------------------------------------- Write Input Parameters. -------------------------------------- !
! ======================================================================================================= !
      
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&neoclassical_input'
         write (unit, '(A, A, A)') '  neo_option = "', trim(neo_option),'"'
         write (unit, '(A, L0)') '  include_neoclassical_terms = ', include_neoclassical_terms
         write (unit, '(A, I0)') '  nradii = ', nradii
         write (unit, '(A, ES0.4)') '  drho = ', drho
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_neoclassical_input

end module namelist_neoclassical_input
