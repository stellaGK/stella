!###############################################################################
!###################### READ PARAMETERS FOR KXKY BOX GRID ######################
!###############################################################################
! Namelist: &parameters_kxky_grids_box
! These flags will allow you to toggle the algorithm choices in stella.
!###############################################################################

module parameters_kxky_grids_box

   implicit none
  
   public :: read_kxky_grids_box
   
   ! Make the namelist public
   public :: read_default_box
   public :: nx, ny, jtwist, jtwistfac, x0, y0, centered_in_rho
   public :: periodic_variation, randomize_phase_shift, phase_shift_angle

   private

   logical :: initialised

   ! Namelist parameters 
   integer :: nx, ny, nalpha, jtwist
   real :: x0, y0, jtwistfac, phase_shift_angle 
   logical :: centered_in_rho, periodic_variation, randomize_phase_shift
   
   ! Parameters set in this module
   integer :: ikx_max, nakx, naky, naky_all
   logical :: reality

contains

   !**********************************************************************
   !                        SET DEFAULT PARAMETERS                       !
   !**********************************************************************
   ! If not specified in the input file these are the default options that 
   ! will be set for all parameters in the namelist &parameters_kxky_grids_box.
   !**********************************************************************
   subroutine read_default_box

      implicit none 

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

   end subroutine read_default_box
    
   !======================================================================
   !================= READ PARAMETERS FOR KXKY BOX GRID ==================
   !======================================================================
   subroutine read_kxky_grids_box(nx_out, ny_out, ikx_max_out, naky_all_out, naky_out, nakx_out, &
         nalpha_out, x0_out, y0_out, jtwist_out, jtwistfac_out, phase_shift_angle_out, &
         centered_in_rho_out, randomize_phase_shift_out, periodic_variation_out, reality_out)

      use mp, only: mp_abort
      use text_options, only: text_option, get_option_value
      use file_utils, only: input_unit, error_unit, input_unit_exist

      implicit none
         
      integer, intent (out) :: nx_out, ny_out, nalpha_out, naky_out, nakx_out
      integer, intent (out) :: jtwist_out, ikx_max_out, naky_all_out
      real, intent (out) :: jtwistfac_out, phase_shift_angle_out, x0_out, y0_out
      logical, intent (out) :: centered_in_rho_out, periodic_variation_out
      logical, intent (out) :: randomize_phase_shift_out, reality_out

      if (initialised) return

      call read_default_box
      call read_input_file_box
      
      ! TODO TODO-HT TODO-GA: Do this in a cleaner way
      nx_out = nx; ny_out = ny; ikx_max_out = ikx_max
      naky_all_out = naky_all; naky_out = naky; nakx_out = nakx
      nalpha_out = nalpha; x0_out = x0; y0_out = y0; jtwist_out = jtwist
      jtwistfac_out = jtwistfac; phase_shift_angle_out = phase_shift_angle
      centered_in_rho_out = centered_in_rho; randomize_phase_shift_out = randomize_phase_shift
      periodic_variation_out = periodic_variation; reality_out = reality

      initialised = .true.
    
   contains

      !**********************************************************************
      !                       READ GRID OPTION FOR KXK                      !
      !**********************************************************************
      ! Read which option to select for the kxky grid layouts
      !**********************************************************************
      subroutine read_input_file_box

         use file_utils, only: input_unit_exist
         use parameters_physics, only: full_flux_surface
         implicit none

         integer :: in_file
         logical :: exist

         namelist /kxky_grids_box/ nx, ny, jtwist, jtwistfac, x0, y0, &
            centered_in_rho, periodic_variation, randomize_phase_shift, phase_shift_angle

         ! note that jtwist and y0 will possibly be modified
         ! later in init_kt_grids_box if they make it out
         ! of this subroutine with negative values
         ! it is necessary to wait until later to do this check
         ! because the values to which they may be set will
         ! depend on information from the geometry module,
         ! which itself may rely on ny from here (number of alphas)

         in_file = input_unit_exist("kxky_grids_box", exist)
         if (exist) read (in_file, nml=kxky_grids_box)
         
         !> Get the number of de-aliased modes in y and x, using reality to halve the number of ky modes
         naky = (ny - 1) / 3 + 1
         nakx = 2 * ((nx - 1) / 3) + 1

         if (full_flux_surface) nalpha = ny
         
         !> Get the ikx index corresponding to kx_max 
         ikx_max = nakx / 2 + 1

         reality = .true.
         !> Get the total number of ky values, including negative ky; 
         !> this is approximately 2/3 ny because ny includes padding to avoid aliasing
         naky_all = 2 * naky - 1
                     
      end subroutine read_input_file_box

   end subroutine read_kxky_grids_box

end module parameters_kxky_grids_box
