!###############################################################################
!################# READ STELLA NAMELISTS FOR THE (KX,KY) GRIDS #################
!###############################################################################
! 
! This module will read the namelists associated with the (kx,ky) grid:
! 
!   kxky_grid_option
!     grid_option = 'default'
!   
!   kxky_grid_range
!     naky = 1.0
!     nakx = 1.0
!     aky_min = 0.0
!     aky_max = 0.0
!     theta0_min = 0.0
!     theta0_max = -1.0
!     akx_min = 0.0
!     akx_max = -1.0
!     kyspacing_option = 'default'
!   
!   kxky_grid_box
!     nx = 1.0
!     ny = 1.0
!     jtwist = -1.0
!     jtwistfac = 1.0
!     x0 = -1.0
!     y0 = -1.0
!     centered_in_rho = .true.
!     periodic_variation = .false.
!     randomize_phase_shift = .false.
!     phase_shift_angle = 0.0
! 
! Text options for <grid_option>:
!    - range: {default, range}
!    - box: {box, annulus, nonlinear}
! 
! Text options for <kyspacing_option>:
!    - linear: {default, linear}
!    - exponential: {exponential}
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
module namelist_kxky_grid

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_kxky_grid_option
   public :: read_namelist_kxky_grid_box
   public :: read_namelist_kxky_grid_range
   
   ! Make parameters accesible to grids_kxky
   public :: kyspacing_linear, kyspacing_exponential
   public :: grid_option_range, grid_option_box
   
   private

   ! Create parameters for <kyspacing_option>
   integer, parameter :: kyspacing_linear = 1
   integer, parameter :: kyspacing_exponential = 2
   
   ! Create parameters for <grid_option>
   integer, parameter :: grid_option_range = 1
   integer, parameter :: grid_option_box = 2
   
   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                           KXKY grid - GRID OPTION                         !
   !****************************************************************************
   subroutine read_namelist_kxky_grid_option(grid_option_switch)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: grid_option_switch

      ! Local variable to set <grid_option_switch>
      character(20) :: grid_option
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call read_default_kxky_grid_option
      call read_input_file_kxky_grid_option

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine read_default_kxky_grid_option

         implicit none

         grid_option = 'default'

      end subroutine read_default_kxky_grid_option

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_kxky_grid_option 

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr

         ! Link text options for <grid_option> to an integer value
         type(text_option), dimension(5), parameter :: grid_options = &
             (/text_option('default', grid_option_range), &
             text_option('range', grid_option_range), &
             text_option('box', grid_option_box), &
             text_option('annulus', grid_option_box), &
             text_option('nonlinear', grid_option_box)/)

         ! Variables in the <kxky_grid_option> namelist
         namelist /kxky_grid_option/ grid_option
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('kxky_grid_option', dexist)
         if (dexist) read (unit=in_file, nml=kxky_grid_option)
         
         ! Read the text option in <grid_option> and store it in <grid_option_switch>
         ierr = error_unit()
         call get_option_value(grid_option, grid_options, grid_option_switch, &
            ierr, 'grid_option in namelist_kxky_grid', .true.)

      end subroutine read_input_file_kxky_grid_option

   end subroutine read_namelist_kxky_grid_option

   !****************************************************************************
   !                          (KX,KY) GRID: BOX OPTION                         !
   !****************************************************************************
   subroutine read_namelist_kxky_grid_box (nx, ny, ikx_max, naky_all, naky, nakx, nalpha, &
      x0, y0, jtwist, jtwistfac, phase_shift_angle, centered_in_rho, randomize_phase_shift, periodic_variation, reality)

      use mp, only: proc0, mp_abort
      use parameters_physics, only: initialised_parameters_physics

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: nx, ny, nalpha
      integer, intent (out) :: jtwist
      integer, intent (out) :: naky, nakx
      integer, intent (out) :: ikx_max, naky_all
      logical, intent (out) :: centered_in_rho, periodic_variation
      logical, intent (out) :: randomize_phase_shift, reality
      real, intent (out) :: jtwistfac
      real, intent (out) :: phase_shift_angle
      real, intent (out) :: x0, y0

      !-------------------------------------------------------------------------
      
      ! Some (kx,ky) grid flags are based on <full_flux_surface>
      ! Therefore, we need to read the physics parameters first
      if (.not. initialised_parameters_physics) then
         call mp_abort('Initialise physics parameters before reading (kx,ky) grid namelists. Aborting.')
      end if

      ! Set the default input parameters and read the input file
      if (.not. proc0) return
      call read_default_kxky_grid_box
      call read_input_file_kxky_grid_box
      call calculate_kxky_grid_box_parameters

   contains

      !------------------------ Default input parameters -----------------------
      subroutine read_default_kxky_grid_box

         implicit none

         ! Default input parameters for the (kx,ky) grids
         nx = 1
         ny = 1
         jtwist = -1
         jtwistfac = 1.
         phase_shift_angle = 0.
         x0 = -1.0
         y0 = -1.0
         nalpha = 1
         centered_in_rho = .true.
         randomize_phase_shift = .false.
         periodic_variation = .false.
         
         ! The following variables can not be modified through the input file
         reality = .true.
         
         ! The following variables are calculated based on (nx,ny)
         nakx = -1
         naky = -1
         ikx_max = -1
         naky_all = -1

         ! Note that jtwist and y0 will possibly be modified
         ! later in init_kxky_grid_box() if they make it out
         ! of this subroutine with negative values.
         ! It is necessary to wait until later to do this check
         ! because the values to which they may be set will
         ! depend on information from the geometry module,
         ! which itself may rely on ny from here (number of alphas)

      end subroutine read_default_kxky_grid_box

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_kxky_grid_box

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables in the <kxky_grid_box> namelist
         namelist /kxky_grid_box/ nx, ny, jtwist, jtwistfac, x0, y0, &
            centered_in_rho, periodic_variation, randomize_phase_shift, phase_shift_angle
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('kxky_grid_box', dexist)
         if (dexist) read (unit=in_file, nml=kxky_grid_box)

      end subroutine read_input_file_kxky_grid_box

      !-------------------------- Calculate parameters -------------------------
      subroutine calculate_kxky_grid_box_parameters
      
         use parameters_physics, only: full_flux_surface

         implicit none

         ! Get the number of de-aliased Fourier modes in x and y, and use the reality condition
         ! (phi[ikx,iky,iz] = conj(phi)[-ikx,-iky,iz]) to halve the number of the ky modes
         ! The number of Fourier modes is approximately 2/3 of (nx,ny) because they include padding to avoid aliasing
         naky = (ny - 1) / 3 + 1
         nakx = 2 * ((nx - 1) / 3) + 1

         ! For a flux tube simulations, we only consider one field line: nalpha = 1
         ! For a full flux surface simulation, we consider ny field lines: nalpha = ny
         if (full_flux_surface) nalpha = ny

         ! Get the ikx index corresponding to kx_max
         ikx_max = nakx / 2 + 1

         ! Get the total number of ky values, including negative ky.
         naky_all = 2 * naky - 1

      end subroutine calculate_kxky_grid_box_parameters

   end subroutine read_namelist_kxky_grid_box

   !****************************************************************************
   !                         (KX,KY) GRID: RANGE OPTION                        !
   !****************************************************************************
   subroutine read_namelist_kxky_grid_range (nalpha, naky, nakx, aky_min, aky_max, &
      akx_min, akx_max, theta0_min, theta0_max, &
      kyspacing_option_switch, phase_shift_angle, ikx_max, naky_all)

      use mp, only: proc0, mp_abort
      use parameters_physics, only: initialised_parameters_physics
      use parameters_physics, only: full_flux_surface

      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: kyspacing_option_switch
      integer, intent (out) :: nalpha, naky, nakx
      integer, intent (out) :: ikx_max, naky_all
      real, intent (out) :: aky_min, aky_max, akx_min, akx_max
      real, intent (out) :: theta0_min, theta0_max
      real, intent (out) :: phase_shift_angle

      ! Local variable to set <kyspacing_option_switch>
      character(20) :: kyspacing_option
      
      !-------------------------------------------------------------------------
      
      ! Some (kx,ky) grid flags are based on <full_flux_surface>
      ! Therefore, we need to read the physics parameters first
      if (.not. initialised_parameters_physics) then
         call mp_abort('Initialise physics parameters before reading (kx,ky) grid namelists. Aborting.')
      end if
         
      ! For full flux surface simulations, we should not be reading the 'range' option for the (kx,ky) grids
      if (full_flux_surface) then
          call mp_abort('The "range" option for the (kx,ky) grids is not supported for full_flux_surface = True. Aborting')
      end if
      
      ! Set the default input parameters and read the input file
      if (.not. proc0) return
      call read_default_kxky_grid_range
      call read_input_file_kxky_grid_range
      call calculate_kxky_grid_range_parameters

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine read_default_kxky_grid_range

         implicit none

         ! Default input parameters for the (kx,ky) grids
         kyspacing_option = 'default'
         nalpha = 1
         naky = 1
         nakx = 1
         aky_min = 0.0
         aky_max = 0.0
         akx_min = 0.0
         akx_max = -1.0
         theta0_min = 0.0
         theta0_max = -1.0
         phase_shift_angle = 0.0
         
         ! The following variables are calculated based on (nakx, naky)
         ikx_max = -1
         naky_all = -1

         ! Note that jtwist and y0 will possibly be modified
         ! later in init_kxky_grid_box() if they make it out
         ! of this subroutine with negative values.
         ! It is necessary to wait until later to do this check
         ! because the values to which they may be set will
         ! depend on information from the geometry module,
         ! which itself may rely on ny from here (number of alphas)

      end subroutine read_default_kxky_grid_range

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_kxky_grid_range

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr

         ! Link text options for <kyspacing_option> to an integer value
         type(text_option), dimension(3), parameter :: kyspacingopts = &
              (/text_option('default', kyspacing_linear), &
              text_option('linear', kyspacing_linear), &
              text_option('exponential', kyspacing_exponential)/)

         ! Variables in the <kxky_grid_range> namelist
         namelist /kxky_grid_range/ naky, nakx, &
              aky_min, aky_max, theta0_min, theta0_max, akx_min, akx_max, kyspacing_option

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('kxky_grid_range', dexist)
         if (dexist) read (unit=in_file, nml=kxky_grid_range)

         ! Read the text option in <kyspacing_option> and store it in <kyspacing_option_switch>
         ierr = error_unit()
         call get_option_value(kyspacing_option, kyspacingopts, kyspacing_option_switch, &
             ierr, 'kyspacing_option in kxky_grid_range', .true.)
         
      end subroutine read_input_file_kxky_grid_range

      !-------------------------- Calculate parameters -------------------------
      subroutine calculate_kxky_grid_range_parameters

         implicit none

         naky_all = naky
         ikx_max = nakx

      end subroutine calculate_kxky_grid_range_parameters

   end subroutine read_namelist_kxky_grid_range

end module namelist_kxky_grid
