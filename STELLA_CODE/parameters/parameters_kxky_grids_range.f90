!###############################################################################
!##################### READ PARAMETERS FOR KXKY RANGE GRID #####################
!###############################################################################
! Namelist: &parameters_kxky_grids_range
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################

module parameters_kxky_grids_range

   implicit none

   public :: read_kxky_grids_range
   
   ! Make the namelist public
   public :: read_default_range
   public :: naky, nakx, aky_min, aky_max, theta0_min, theta0_max
   public :: akx_min, akx_max, kyspacing_option

   private
   
   ! Variables in the namelist
   integer :: naky, nakx
   real :: aky_min, aky_max, akx_min, akx_max
   real :: theta0_min, theta0_max
   character(20) :: kyspacing_option
   
   ! Parameters set in this module
   integer :: ikx_max, kyspacing_option_switch, naky_all, nalpha
   real :: phase_shift_angle

   ! Local variable
   logical :: initialised

contains

   !**********************************************************************
   !                        SET DEFAULT PARAMETERS                       !
   !**********************************************************************
   ! If not specified in the input file these are the default options that 
   ! will be set for all parameters under the namelist 
   ! &numerical'.
   !**********************************************************************
   subroutine read_default_range 

      implicit none

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

   !======================================================================
   !================= READ PARAMETERS FOR KXKY range GRID ================
   !======================================================================
   subroutine read_kxky_grids_range (nalpha_out, naky_out, nakx_out, aky_min_out, aky_max_out, &
           akx_min_out, akx_max_out, theta0_min_out, theta0_max_out, &
           kyspacing_option_switch_out, phase_shift_angle_out, ikx_max_out, naky_all_out)

      use mp, only: mp_abort

      implicit none

      integer, intent (out) :: nalpha_out, naky_out, nakx_out
      real, intent (out) :: aky_min_out, aky_max_out, akx_min_out, akx_max_out
      real, intent (out) :: theta0_min_out, theta0_max_out
      integer, intent (out) :: kyspacing_option_switch_out
      real, intent (out) :: phase_shift_angle_out
      integer, intent (out) :: ikx_max_out, naky_all_out

      integer, parameter :: kyspacing_linear = 1, kyspacing_exponential = 2

      if (initialised) return

      call read_default_range 
      call read_input_file_range 
      
      ! TODO TODO-HT TODO-GA: Do this in a cleaner way
      nalpha_out = nalpha; naky_out = naky; nakx_out = nakx
      aky_min_out = aky_min; aky_max_out = aky_max
      akx_min_out = akx_min; akx_max_out = akx_max
      theta0_min_out = theta0_min; theta0_max_out = theta0_max
      kyspacing_option_switch_out = kyspacing_option_switch
      phase_shift_angle_out = phase_shift_angle
      ikx_max_out = ikx_max; naky_all_out = naky_all

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
