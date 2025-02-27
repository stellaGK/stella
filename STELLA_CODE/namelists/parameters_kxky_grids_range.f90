!###############################################################################
!##################### READ PARAMETERS FOR KXKY RANGE GRID #####################
!###############################################################################
! Namelist: &parameters_kxky_grids_range
! These flags will allow you to toggle the algorithm choices in stella.
! 
! Note that we do not want to make the variables from the <kxky_grids_box>
! namelist public, because only the variables from parameters_kxky_grids.f90 
! should be used throughout stella, this script only handles the namelist.
!###############################################################################

module parameters_kxky_grids_range

   implicit none

   ! parameters_kxky_grids.f90 needs read_kxky_grids_range
   ! update_input_file.f90 needs read_default_range
   public :: read_kxky_grids_range 
   public :: read_default_range

   private

   logical :: initialised

contains

   !============================================================================
   !========================= SET DEFAULT PARAMETERS ==========================!
   !============================================================================
   ! If not specified in the input file these are the default options that 
   ! will be set for all parameters under the namelist &kxky_grids_range.
   !============================================================================
   subroutine read_default_range(naky, nakx, aky_min, aky_max, &
         akx_min, akx_max, theta0_min, theta0_max, kyspacing_option)

      implicit none

      integer, intent (out) :: naky, nakx
      real, intent (out) :: aky_min, aky_max, akx_min, akx_max
      real, intent (out) :: theta0_min, theta0_max
      character(20), intent (out) :: kyspacing_option

      naky = 1
      nakx = 1
      aky_min = 0.0
      aky_max = 0.0
      
      !> set these to be nonsense values
      !> so can check later if they've been set
      akx_min = 0.0
      akx_max = -1.0
      theta0_min = 0.0
      theta0_max = -1.0
      kyspacing_option = 'default'

   end subroutine read_default_range

   !============================================================================
   !==================== READ PARAMETERS FOR KXKY range GRID ===================
   !============================================================================
   subroutine read_kxky_grids_range(nalpha, naky, nakx, &
           aky_min, aky_max, akx_min, akx_max, theta0_min, theta0_max, &
           kyspacing_option_switch, phase_shift_angle, ikx_max, naky_all)

      use mp, only: mp_abort

      implicit none

      integer, intent (out) :: nalpha, naky, nakx
      real, intent (out) :: aky_min, aky_max, akx_min, akx_max
      real, intent (out) :: theta0_min, theta0_max
      integer, intent (out) :: kyspacing_option_switch
      real, intent (out) :: phase_shift_angle
      integer, intent (out) :: ikx_max, naky_all

      ! Local variables
      integer, parameter :: kyspacing_linear = 1, kyspacing_exponential = 2
      character(20) :: kyspacing_option

      if (initialised) return

      call read_default_range(naky, nakx, aky_min, aky_max, &
         akx_min, akx_max, theta0_min, theta0_max, kyspacing_option)
      call read_input_file_range 

      initialised = .true.

   contains

      !**********************************************************************
      !                   READ INPUT FOR KXKY GRIDS, RANGE                  !
      !**********************************************************************
      ! Read which option to select for the kxky grid layouts
      !**********************************************************************
      subroutine read_input_file_range

         use file_utils, only: input_unit, error_unit, input_unit_exist
         use parameters_physics, only: full_flux_surface
         use text_options, only: text_option, get_option_value

         implicit none

         type(text_option), dimension(3), parameter :: kyspacingopts = &
              (/text_option('default', kyspacing_linear), &
              text_option('linear', kyspacing_linear), &
              text_option('exponential', kyspacing_exponential)/) 

         integer :: ierr, in_file
         logical :: exist

         namelist /kxky_grids_range/ naky, nakx, aky_min, aky_max, &
              theta0_min, theta0_max, akx_min, akx_max, kyspacing_option
      
         ! Phase shift angle and nalpha are not part of the namelist!
         phase_shift_angle = 0.
         nalpha = 1

         ! note that jtwist and y0 will possibly be modified
         ! later in init_kt_grids_range if they make it out
         ! of this subroutine with negative values
         ! it is necessary to wait until later to do this check
         ! because the values to which they may be set will
         ! depend on information from the geometry module,
         ! which itself may rely on ny from here (number of alphas)

         in_file = input_unit_exist("kxky_grids_range", exist)
         if (exist) read (in_file, nml=kxky_grids_range)
         
         if (full_flux_surface) then
             write (*, *) '!!! ERROR !!!'
             write (*, *) 'kt_grids "range" option is not supported for full_flux_surface = T. aborting'
             write (*, *) '!!! ERROR !!!'
             call mp_abort('kt_grids "range" option is not supported for full_flux_surface = T. aborting')
         end if
         
         ierr = error_unit()
         call get_option_value(kyspacing_option, kyspacingopts, kyspacing_option_switch, &
             ierr, "kyspacing_option in kt_grids_range_parameters", .true.)

         naky_all = naky
         ikx_max = nakx
         
      end subroutine read_input_file_range

    end subroutine read_kxky_grids_range

end module parameters_kxky_grids_range
