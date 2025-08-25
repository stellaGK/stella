!###############################################################################
!############ READ STELLA NAMELISTS FOR RADIAL VARIATION (MULTIBOX) ############
!###############################################################################
! 
! This module will read the namelists associated with radial variation:
! 
!   multibox_parameters
!     ky_solve_radial = 0.0
!     ky_solve_real = .false.
!     include_geometric_variation = .true.
!     include_pressure_variation = .false.
!     smooth_zf = .false.
!     zf_option = 'default'
!     lr_debug_option = 'default'
!     krook_option = 'default'
!     krook_size = 0.0
!     rk_step = .false.
!     nu_krook_mb = 0.0
!     mb_debug_step = -1.0
!     krook_exponent = 0.0
!     krook_efold = 3.0
!     phi_bound = 0.0
!     phi_pow = 0.0
!     use_dirichlet_bc = .false.
!     boundary_size = 4.0
!     krook_size = 0.0
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
module namelist_radial_variation

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_radial_variation

   ! Parameters need to be public (krook_option)
   public :: krook_option_default, krook_option_flat
   public :: krook_option_linear, krook_option_exp, krook_option_exp_rev
   
   ! Parameters need to be public (zf_option)
   public :: mb_zf_option_default, mb_zf_option_skip_ky0
   public ::  mb_zf_option_zero_ky0, mb_zf_option_zero_fsa
         
   ! Parameters need to be public (lr_debug_option)
   public :: lr_debug_option_default, lr_debug_option_L, lr_debug_option_R
         
   private

   ! Create parameters for <krook_option>
   integer, parameter :: krook_option_default = 2
   integer, parameter :: krook_option_flat = 0
   integer, parameter :: krook_option_linear = 1
   integer, parameter :: krook_option_exp = 2
   integer, parameter :: krook_option_exp_rev = 3
   
   ! Create parameters for <zf_option>
   integer, parameter :: mb_zf_option_default = 0
   integer, parameter :: mb_zf_option_skip_ky0 = 1
   integer, parameter :: mb_zf_option_zero_ky0 = 2
   integer, parameter :: mb_zf_option_zero_fsa = 3
   
   ! Create parameters for <lr_debug_option>
   integer, parameter:: lr_debug_option_default = 0
   integer, parameter:: lr_debug_option_L = 1
   integer, parameter:: lr_debug_option_R = 2

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                               RADIAL VARIATION                            !
   !****************************************************************************
   subroutine read_namelist_radial_variation(ky_solve_real, ky_solve_radial, &
      include_pressure_variation, include_geometric_variation, &
      smooth_zf, lr_debug_switch, krook_option_switch, mb_zf_option_switch, &
      rk_step, nu_krook_mb, mb_debug_step, &
      krook_exponent, krook_efold, phi_bound, phi_pow, &
      use_dirichlet_bc, boundary_size, krook_size)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      logical, intent (out) :: ky_solve_real
      logical, intent (out) :: include_pressure_variation
      logical, intent (out) :: include_geometric_variation
      logical, intent (out) :: smooth_zf
      logical, intent (out) :: rk_step
      logical, intent (out) :: use_dirichlet_bc
      real, intent (out) ::  nu_krook_mb
      real, intent (out) :: krook_exponent, krook_efold
      real, intent (out) :: phi_bound, phi_pow
      integer, intent (out) :: ky_solve_radial
      integer, intent (out) :: mb_debug_step
      integer, intent (out) :: boundary_size, krook_size
      integer, intent (out) :: lr_debug_switch, krook_option_switch, mb_zf_option_switch

      ! Local variables to set <krook_option_switch>, <mb_zf_option_switch> and <lr_debug_switch>
      character(30) :: zf_option, krook_option, lr_debug_option
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_radial_variation
      call read_input_file_radial_variation
      call check_namelist_radial_variation

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_radial_variation

         implicit none

         ky_solve_radial = 0
         ky_solve_real =.false.
         include_pressure_variation = .false.
         include_geometric_variation = .false.
         smooth_zf = .false.
         lr_debug_option = 'default'
         krook_option = 'default'
         zf_option = 'default'
         rk_step = .false.
         nu_krook_mb = 0.0
         mb_debug_step = -1.0
         krook_exponent = 0.0
         krook_efold = 3.0
         phi_bound = 0.0
         phi_pow = 0.0
         use_dirichlet_bc = .false.
         boundary_size = 4.0
         krook_size = 0.0
         
      end subroutine set_default_parameters_radial_variation

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_radial_variation

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr
         
         ! Link text options for <krook_option> to an integer value
         type(text_option), dimension(5), parameter :: krook_opts = &
            (/text_option('default', krook_option_default), &
              text_option('flat', krook_option_flat), &
              text_option('linear', krook_option_linear), &
              text_option('exp', krook_option_exp), &
              text_option('exp_reverse', krook_option_exp_rev)/)
         
         ! Link text options for <mb_zf_option> to an integer value
         type(text_option), dimension(4), parameter :: mb_zf_opts = &
            (/text_option('default', mb_zf_option_default), &
              text_option('skip_ky0', mb_zf_option_skip_ky0), &
              text_option('zero_ky0', mb_zf_option_zero_ky0), &
              text_option('zero_fsa', mb_zf_option_zero_fsa)/)
         
         ! Link text options for <lr_debug_option> to an integer value
         type(text_option), dimension(3), parameter :: lr_db_opts = &
            (/text_option('default', lr_debug_option_default), &
              text_option('L', lr_debug_option_L), &
              text_option('R', lr_debug_option_R)/)

         ! Variables in the <multibox_parameters> namelist
         namelist /multibox_parameters/ ky_solve_real, ky_solve_radial, &
            include_pressure_variation, include_geometric_variation, &
            smooth_zf, lr_debug_option, krook_option, zf_option, &
            rk_step, nu_krook_mb, mb_debug_step, &
            krook_exponent, krook_efold, phi_bound, phi_pow, &
            use_dirichlet_bc, boundary_size, krook_size
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('multibox_parameters', dexist)
         if (dexist) read (unit=in_file, nml=multibox_parameters)

         ! Read the text option and store it the corresponding switch
         ierr = error_unit()
         call get_option_value(krook_option, krook_opts, krook_option_switch, &
            ierr, 'krook_option in multibox_parameters')
         call get_option_value(zf_option, mb_zf_opts, mb_zf_option_switch, &
            ierr, 'zf_option in multibox_parameters')
         call get_option_value(lr_debug_option, lr_db_opts, lr_debug_switch, &
            ierr, 'lr_debug_option in multibox_parameters')

      end subroutine read_input_file_radial_variation

      !------------------------- Check input parameters ------------------------
      subroutine check_namelist_radial_variation

         implicit none
         
         if (krook_size > boundary_size) krook_size = boundary_size

      end subroutine check_namelist_radial_variation

   end subroutine read_namelist_radial_variation

end module namelist_radial_variation
