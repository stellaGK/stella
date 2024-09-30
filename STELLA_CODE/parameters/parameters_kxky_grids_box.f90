!###############################################################################
!###################### READ PARAMETERS FOR KXKY BOX GRID ######################
!###############################################################################
! Namelist: &parameters_kxky_grids_box
! These flags will allow you to toggle the algorithm choices in stella.
! 
! Note that we do not want to make the variables from the <kxky_grids_box>
! namelist public, because only the variables from parameters_kxky_grids.f90 
! should be used throughout stella, this script only handles the namelist.
!###############################################################################

module parameters_kxky_grids_box

   implicit none
  
   ! parameters_kxky_grids.f90 needs read_kxky_grids_box
   ! update_input_file.f90 needs read_default_box
   public :: read_kxky_grids_box 
   public :: read_default_box

   private

   logical :: initialised

contains

   !============================================================================
   !========================= SET DEFAULT PARAMETERS ==========================!
   !============================================================================
   ! If not specified in the input file these are the default options that 
   ! will be set for all parameters in the namelist &parameters_kxky_grids_box.
   !============================================================================
   subroutine read_default_box(nx, ny, nalpha, &
         jtwist, x0, y0, jtwistfac, phase_shift_angle, &
         centered_in_rho, periodic_variation, randomize_phase_shift)

      implicit none 
      
      integer, intent (out) :: nx, ny, nalpha, jtwist
      real, intent (out) :: jtwistfac, phase_shift_angle, x0, y0
      logical, intent (out) :: centered_in_rho, periodic_variation
      logical, intent (out) :: randomize_phase_shift

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
    
   !============================================================================
   !==================== READ PARAMETERS FOR KXKY BOX GRID =====================
   !============================================================================
   subroutine read_kxky_grids_box(nx, ny, ikx_max, naky_all, naky, nakx, &
         nalpha, x0, y0, jtwist, jtwistfac, phase_shift_angle, &
         centered_in_rho, randomize_phase_shift, periodic_variation, reality)

      use mp, only: mp_abort
      use text_options, only: text_option, get_option_value
      use file_utils, only: input_unit, error_unit, input_unit_exist

      implicit none
         
      integer, intent (out) :: nx, ny, nalpha, naky, nakx
      integer, intent (out) :: jtwist, ikx_max, naky_all
      real, intent (out) :: jtwistfac, phase_shift_angle, x0, y0
      logical, intent (out) :: centered_in_rho, periodic_variation
      logical, intent (out) :: randomize_phase_shift, reality

      if (initialised) return

      call read_default_box(nx, ny, nalpha, &
         jtwist, x0, y0, jtwistfac, phase_shift_angle, &
         centered_in_rho, periodic_variation, randomize_phase_shift)
      call read_input_file_box()

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
