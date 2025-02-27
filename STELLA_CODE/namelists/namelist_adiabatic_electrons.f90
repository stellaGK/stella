!###############################################################################
!###################### READ ADIABATIC ELECTRONS NAMELIST ######################
!###############################################################################
! Read the <adiabatic_electrons> namelist in the input file:
! 
! &adiabatic_electroms
!  adiabatic_option = 'field-line-average-term'
!  tite = 1.0
!  nine = 1.0
!/
!###############################################################################
module namelist_adiabatic_electrons

   implicit none

   public :: read_namelist
   
   private

contains

   subroutine read_namelist(adiabatic_option_switch, nine, tite, &
      adiabatic_option_fieldlineavg, adiabatic_option_periodic)

      use mp, only: proc0

      implicit none

      real, intent(out) :: tite, nine 
      integer, intent(out) :: adiabatic_option_switch
      integer, intent(out) :: adiabatic_option_periodic
      integer, intent(out) :: adiabatic_option_fieldlineavg
      
      ! Local variable to set <adiabatic_option_switch>
      character(30) :: adiabatic_option

      if (.not. proc0) return
      call set_default_parameters
      call read_input_file
      call check_inputs

   contains

      !**********************************************************************
      !                         DEFAULT PARAMETERS                          !
      !**********************************************************************
      ! Set the default input parameters.
      !**********************************************************************
      subroutine set_default_parameters()

         implicit none

         adiabatic_option = 'field-line-average-term'
         tite = 1.0
         nine = 1.0

      end subroutine set_default_parameters


      !**********************************************************************
      !                          READ INPUT FILE                            !
      !**********************************************************************
      ! Overwrite any default options with those specified in the input file. 
      ! Then change the other parameters consistently.
      !**********************************************************************
      subroutine read_input_file

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr, in_file
         logical :: dexist
         
         ! Create parameters for the text option
         integer, parameter :: periodic = 1 
         integer, parameter :: fieldlineavg = 2
         
         ! Link text options for <adiabatic_option> to an integer value
         type(text_option), dimension(6), parameter :: adiabaticopts = &
            (/text_option('default', fieldlineavg), &
            text_option('no-field-line-average-term', periodic), &
            text_option('field-line-average-term', fieldlineavg), &
            text_option('iphi00=0', periodic), &
            text_option('iphi00=1', periodic), &
            text_option('iphi00=2', fieldlineavg)/)

         ! Define variables in the <adiabatic_electrons> namelist
         namelist /adiabatic_electrons/ adiabatic_option, tite, nine
         
         !----------------------------------------------------------------------
         
         ! Save the values of the adiabatic option
         adiabatic_option_periodic = periodic
         adiabatic_option_fieldlineavg = fieldlineavg

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist("adiabatic_electrons", dexist)
         if (dexist) read (unit=in_file, nml=adiabatic_electrons)
         
         ! Read the text option in <adiabatic_option> and store it in <adiabatic_option_switch>
         ierr = error_unit()
         call get_option_value(adiabatic_option, adiabaticopts, adiabatic_option_switch, &
            ierr, "adiabatic_option in parameters_physics")

      end subroutine read_input_file


      !**********************************************************************
      !                           CHECK INPUTS                             !
      !**********************************************************************
      ! Make sure that the input variables are set correctly.
      !**********************************************************************
      subroutine check_inputs

         implicit none

      end subroutine check_inputs

   end subroutine read_namelist

end module namelist_adiabatic_electrons

