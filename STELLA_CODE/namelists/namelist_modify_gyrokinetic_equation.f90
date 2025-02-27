!###############################################################################
!################## READ MODIFY GYROKINETIC EQUATION NAMELIST ##################
!###############################################################################
! Read the <modify_gyrokinetic_equation> namelist in the input file:
! 
! &modify_gyrokinetic_equation
!  suppress_zonal_interaction = .false.
!/
!###############################################################################
module namelist_modify_gyrokinetic_equation

   implicit none

   public :: read_namelist
   
   private

contains

   subroutine read_namelist(suppress_zonal_interaction)

      use mp, only: proc0

      implicit none

      logical, intent(out) :: suppress_zonal_interaction

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

         suppress_zonal_interaction = .false.

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
         integer :: in_file
         logical :: dexist 

         ! Define variables in the <modify_gyrokinetic_equation> namelist
         namelist /modify_gyrokinetic_equation/ suppress_zonal_interaction

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist("modify_gyrokinetic_equation", dexist)
         if (dexist) read (unit=in_file, nml=modify_gyrokinetic_equation)

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

end module namelist_modify_gyrokinetic_equation

