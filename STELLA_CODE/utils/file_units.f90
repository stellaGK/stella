!###############################################################################
!################################## FILE UNITS #################################
!###############################################################################
! Save the unit of the error file so we don't need to call a function each
! time we need it. And save the unit number of the default input file, since
! many modules in the "read_namelists" folder will write to it.
! 
! Note that the fixed unit numbers should really be using the 
! call get_unused_unit(geofile_unit) function instead.
!###############################################################################
module file_units

   implicit none
   
   ! Make the routines available to other modules
   public :: init_file_units
   
   ! Make the routines available to other modules
   public :: open_input_file_with_defaults
   public :: close_input_file_with_defaults

   ! Make unit numbers available to other routines
   public :: unit_error_file
   public :: unit_sfincs_file
   public :: unit_input_file_with_defaults
   public :: unit_euterpe_input
   public :: unit_euterpe_output
   public :: unit_response_matrix

   private

   ! Unit numbers for opening files
   integer :: unit_error_file
   integer :: unit_input_file_with_defaults
   
   ! Some units just have fixed values
   integer, parameter :: unit_sfincs_file = 999
   integer, parameter :: unit_euterpe_input = 1099
   integer, parameter :: unit_euterpe_output = 1098
   integer, parameter :: unit_response_matrix = 70

contains

   !****************************************************************************
   !                           INIT FIXED FILE UNITS                           !
   !****************************************************************************
   subroutine init_file_units
      
      use file_utils, only: error_unit
   
      implicit none
      
      !-------------------------------------------------------------------------
      
      ! Get the unit of the error file
      unit_error_file = error_unit()
      
   end subroutine init_file_units
   
   
   !****************************************************************************
   !                            OPEN AND CLOSE FILES                           !
   !****************************************************************************
   ! Write the input parameters from the input files, and the default input
   ! parameters used in stella, to *.input.
   !****************************************************************************
      
   !---------------------- Open input file with defaults -----------------------
   subroutine open_input_file_with_defaults
   
      use file_utils, only: open_output_file
      
      call open_output_file(unit_input_file_with_defaults, '.input')
   
   end subroutine open_input_file_with_defaults

   !---------------------- Close input file with defaults -----------------------
   subroutine close_input_file_with_defaults

      implicit none

      close (unit_input_file_with_defaults)

   end subroutine close_input_file_with_defaults
   
end module file_units
