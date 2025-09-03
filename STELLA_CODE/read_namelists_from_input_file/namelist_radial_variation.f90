!###############################################################################
!############ READ STELLA NAMELISTS FOR RADIAL VARIATION (MULTIBOX) ############
!###############################################################################
! 
! Note that the sources are sinks, and flow shear, have been implemented for the
! radially global version of stella, discussed in  [2022 - St-Onge - A novel approach
! to radially global gyrokinetic simulation using the flux-tube code stella]
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
!   sources
!     source_option = 'none'
!     nu_krook = 0.05
!     tcorr_source = 0.02
!     ikxmax_source = 1.0
!     krook_odd = .true.
!     exclude_boundary_regions = .false.
!     tcorr_source_qn = 0.0
!     exclude_boundary_regions_qn = .false.
!     from_zero = .true.
!     conserve_momentum = .false.
!     conserve_density = .false.
!   
!   flow_shear
!     prp_shear_enabled = .false.
!     hammett_flow_shear = .true.
!     g_exb = 0.0
!     g_exbfac = 1.0
!     omprimfac = 1.0
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
   public :: read_namelist_multibox
   public :: read_namelist_sources
   public :: read_namelist_flow_shear

   ! Parameters need to be public (krook_option)
   public :: krook_option_default, krook_option_flat
   public :: krook_option_linear, krook_option_exp, krook_option_exp_rev
   
   ! Parameters need to be public (zf_option)
   public :: mb_zf_option_default, mb_zf_option_skip_ky0
   public ::  mb_zf_option_zero_ky0, mb_zf_option_zero_fsa
         
   ! Parameters need to be public (lr_debug_option)
   public :: lr_debug_option_default, lr_debug_option_L, lr_debug_option_R

   ! Parameters need to be public (source_option)
   public :: source_option_none, source_option_krook, source_option_projection
         
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
   
   ! Create parameters for <source_option>
   integer, parameter :: source_option_none = 1
   integer, parameter :: source_option_krook = 2
   integer, parameter :: source_option_projection = 3

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                                  MULTI BOX                                !
   !****************************************************************************
   subroutine read_namelist_multibox(ky_solve_real, ky_solve_radial, &
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
      call set_default_parameters_multibox
      call read_input_file_multibox
      call check_namelist_multibox
      call write_parameters_to_input_file_multibox

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_multibox

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
         
      end subroutine set_default_parameters_multibox

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_multibox

         use file_units, only: unit_error_file
         use file_utils, only: input_unit_exist
         use text_options, only: text_option, get_option_value

         implicit none
         
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
         call get_option_value(krook_option, krook_opts, krook_option_switch, &
            unit_error_file, 'krook_option in multibox_parameters')
         call get_option_value(zf_option, mb_zf_opts, mb_zf_option_switch, &
            unit_error_file, 'zf_option in multibox_parameters')
         call get_option_value(lr_debug_option, lr_db_opts, lr_debug_switch, &
            unit_error_file, 'lr_debug_option in multibox_parameters')

      end subroutine read_input_file_multibox

      !------------------------- Check input parameters ------------------------
      subroutine check_namelist_multibox

         implicit none
         
         if (krook_size > boundary_size) krook_size = boundary_size

      end subroutine check_namelist_multibox
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file_multibox

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&multibox_parameters'
         write (unit, '(A, A, A)') '  lr_debug_option = "', lr_debug_option, '"'
         write (unit, '(A, A, A)') '  krook_option = "', krook_option, '"'
         write (unit, '(A, A, A)') '  zf_option = "', zf_option, '"'
         write (unit, '(A, L0)') '  ky_solve_real = ', ky_solve_real
         write (unit, '(A, I0)') '  ky_solve_radial = ', ky_solve_radial
         write (unit, '(A, L0)') '  include_pressure_variation = ', include_pressure_variation
         write (unit, '(A, L0)') '  include_geometric_variation = ', include_geometric_variation
         write (unit, '(A, L0)') '  smooth_zf = ', smooth_zf
         write (unit, '(A, L0)') '  rk_step = ', rk_step
         write (unit, '(A, ES0.4)') '  nu_krook_mb = ', nu_krook_mb
         write (unit, '(A, I0)') '  mb_debug_step = ', mb_debug_step
         write (unit, '(A, ES0.4)') '  krook_exponent = ', krook_exponent
         write (unit, '(A, ES0.4)') '  krook_efold = ', krook_efold
         write (unit, '(A, ES0.4)') '  phi_bound = ', phi_bound
         write (unit, '(A, ES0.4)') '  phi_pow = ', phi_pow
         write (unit, '(A, L0)') '  use_dirichlet_bc = ', use_dirichlet_bc
         write (unit, '(A, I0)') '  boundary_size = ', boundary_size
         write (unit, '(A, I0)') '  krook_size = ', krook_size
         write (unit, '(A)') '/'
         write (unit, '(A)') ''
      
      end subroutine write_parameters_to_input_file_multibox

   end subroutine read_namelist_multibox

   !****************************************************************************
   !                                  SOURCES                                  !
   !****************************************************************************
   ! Some sources flags are based on <radial_variation> and <periodic_variation>
   ! Therefore, we need to read the physics parameters and (kx,ky) grids first
   !****************************************************************************
   subroutine read_namelist_sources(periodic_variation, source_option_switch, nu_krook, tcorr_source, &
      ikxmax_source, krook_odd, exclude_boundary_regions, &
      tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
      conserve_momentum, conserve_density)

      use mp, only: proc0, mp_abort

      implicit none
      
      ! Input parameters from other namelists
      logical, intent (out) :: periodic_variation

      ! Variables that are read from the input file
      integer, intent (out) :: source_option_switch, ikxmax_source
      logical, intent (out) :: krook_odd, exclude_boundary_regions
      logical, intent (out) :: exclude_boundary_regions_qn, from_zero
      logical, intent (out) :: conserve_momentum, conserve_density
      real, intent (out) :: nu_krook, tcorr_source
      real, intent (out) :: tcorr_source_qn
      
      ! Local variable to set <source_option>
      character(30) :: source_option
      
      !-------------------------------------------------------------------------

      ! Set the default input parameters and read the input file
      if (.not. proc0) return
      call set_default_parameters_sources
      call read_input_file_sources
      call check_inputs_sources
      call write_parameters_to_input_file_sources

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_sources

         use parameters_physics, only: radial_variation

         implicit none

         ! Default parameters for the sources
         source_option = 'none'
         nu_krook = 0.05
         tcorr_source = 0.02
         tcorr_source_qn = 0.0
         from_zero = .true.
         conserve_momentum = .false.
         conserve_density = .false.
         exclude_boundary_regions = radial_variation .and. .not. periodic_variation
         
         ! Damp only the odd mode that can affect profiles
         krook_odd = .true.
         
         ! For periodic variation we use kx=0 and kx=1
         ikxmax_source = 1.0
         if (periodic_variation) ikxmax_source = 2 

      end subroutine set_default_parameters_sources

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_sources

         use file_utils, only: input_unit_exist
         use file_units, only: unit_error_file
         use text_options, only: text_option, get_option_value

         implicit none
         
         ! Link text options for <source_option> to an integer value
         type(text_option), dimension(4), parameter :: sourceopts = &
            (/text_option('default', source_option_none), &
              text_option('none', source_option_none), &
              text_option('krook', source_option_krook), &
              text_option('projection', source_option_projection)/)

         ! Variables in the <sources> namelist
         namelist /sources/ source_option, nu_krook, tcorr_source, &
            ikxmax_source, krook_odd, exclude_boundary_regions, &
            tcorr_source_qn, exclude_boundary_regions_qn, from_zero, &
            conserve_momentum, conserve_density
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('sources', dexist)
         if (dexist) read (unit=in_file, nml=sources)

         ! Read the text option in <source_option> and store it in <source_option_switch>
         call get_option_value(source_option, sourceopts, source_option_switch, &
            unit_error_file, 'source_option in sources')

      end subroutine read_input_file_sources

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_sources

         implicit none

         if (tcorr_source_qn < 0) tcorr_source_qn = tcorr_source

      end subroutine check_inputs_sources
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file_sources

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------
         
         ! Don't print values if no sources are selected
         if (source_option_switch == source_option_none) return

         ! Print input parameters
         write (unit, '(A)') '&sources'
         write (unit, '(A, A, A)') '  source_option = "', trim(source_option), '"'
         write (unit, '(A, ES0.4)') '  nu_krook = ', nu_krook
         write (unit, '(A, ES0.4)') '  tcorr_source = ', tcorr_source
         write (unit, '(A, I0)') '  ikxmax_source = ', ikxmax_source
         write (unit, '(A, L0)') '  krook_odd = ', krook_odd
         write (unit, '(A, L0)') '  exclude_boundary_regions = ', exclude_boundary_regions
         write (unit, '(A, ES0.4)') '  tcorr_source_qn = ', tcorr_source_qn
         write (unit, '(A, L0)') '  exclude_boundary_regions_qn = ', exclude_boundary_regions_qn
         write (unit, '(A, L0)') '  from_zero = ', from_zero
         write (unit, '(A, L0)') '  conserve_momentum = ', conserve_momentum
         write (unit, '(A, L0)') '  conserve_density = ', conserve_density
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file_sources

   end subroutine read_namelist_sources
   
   !****************************************************************************
   !                               FLOW SHEAR TERMS                            !
   !****************************************************************************
   subroutine read_namelist_flow_shear(prp_shear_enabled, hammett_flow_shear, &
      g_exb, g_exbfac, omprimfac)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      logical, intent(out) :: prp_shear_enabled, hammett_flow_shear 
      real :: g_exb, g_exbfac, omprimfac
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_flow_shear
      call read_input_file_flow_shear
      call write_parameters_to_input_file

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_flow_shear

         implicit none

         prp_shear_enabled = .false.
         hammett_flow_shear = .true.
         g_exb = 0.0
         g_exbfac = 1.0
         omprimfac = 1.0

      end subroutine set_default_parameters_flow_shear

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_flow_shear

         use file_utils, only: input_unit_exist
         implicit none

         namelist /flow_shear/ prp_shear_enabled, hammett_flow_shear, g_exb, g_exbfac, omprimfac
         in_file = input_unit_exist('flow_shear', dexist)
         if (dexist) read (unit=in_file, nml=flow_shear)

      end subroutine read_input_file_flow_shear
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&flow_shear'
         write (unit, '(A, L0)') '  prp_shear_enabled = ', prp_shear_enabled
         write (unit, '(A, L0)') '  hammett_flow_shear = ', hammett_flow_shear
         write (unit, '(A, ES0.4)') '  g_exb = ', g_exb
         write (unit, '(A, ES0.4)') '  g_exbfac = ', g_exbfac
         write (unit, '(A, ES0.4)') '  omprimfac = ', omprimfac
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_flow_shear
   
end module namelist_radial_variation
