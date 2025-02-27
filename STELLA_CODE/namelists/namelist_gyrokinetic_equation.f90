!###############################################################################
!##################### READ GYROKINETIC EQUATION NAMELIST ######################
!###############################################################################
! Read the <gyrokinetic_equation> namelist in the input file:
! 
! &gyrokinetic_equation
!  include_parallel_nonlinearity = .false.
!  include_parallel_streaming = .true.
!  include_mirror = .true.
!  nonlinear = .false.
!  xdriftknob = 1.0
!  ydriftknob = 1.0
!  wstarknob = 1.0
!/
!###############################################################################
module namelist_gyrokinetic_equation

   implicit none

   public :: read_namelist
   
   private

contains

   subroutine read_namelist(xdriftknob, ydriftknob, wstarknob, nonlinear, &
      include_mirror, include_parallel_streaming, include_parallel_nonlinearity)

      use mp, only: proc0

      implicit none

      logical, intent(out) :: include_parallel_nonlinearity
      logical, intent(out) :: include_parallel_streaming
      logical, intent(out) :: include_mirror, nonlinear
      real, intent(out) :: xdriftknob, ydriftknob, wstarknob

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

         include_parallel_nonlinearity = .false.
         include_parallel_streaming = .true.
         include_mirror = .true.
         nonlinear = .false.
         xdriftknob = 1.0
         ydriftknob = 1.0
         wstarknob = 1.0

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

         ! Define variables in the <gyrokinetic_equation> namelist
         namelist /gyrokinetic_equation/ xdriftknob, ydriftknob, wstarknob, nonlinear, &
            include_mirror, include_parallel_streaming, include_parallel_nonlinearity

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist("gyrokinetic_equation", dexist)
         if (dexist) read (unit=in_file, nml=gyrokinetic_equation)

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

end module namelist_gyrokinetic_equation

