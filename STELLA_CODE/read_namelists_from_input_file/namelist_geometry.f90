!###############################################################################
!###################### READ STELLA NAMELISTS FOR GEOMETRY #####################
!###############################################################################
! 
! This module will read the namelists associated with geometry:
! 
!   geometry_options
!     geometry_option = 'local'
!     q_as_x = .false.
!   
!   geometry_miller
!     rhoc = 0.5
!     rmaj = 3.0
!     shift = 0.0
!     qinp = 1.4
!     shat = 0.8
!     kappa = 0.0
!     kapprim = 0.0
!     tri = 0.0
!     triprim = 0.0
!     rgeo = 3.0
!     betaprim = 0.0
!     betadbprim = 0.0
!     d2qdr2 = 0.0
!     d2psidr2 = 0.0
!     nzed_local = 128.0
!     read_profile_variation = .false.
!     write_profile_variation = .false.
!   
!   geometry_vmec
!     alpha0 = 0.0
!     zeta_center = 0.0
!     rectangular_cross_section = .false.
!     nfield_periods = -1.0
!     torflux = 0.6354167
!     z_grid_refinement_factor = 1.0
!     surface_option = 0.0
!     radial_coordinate = 'sgn(psi_t)psi_t'
!     verbose = .true.
!     vmec_filename = 'equilibria/wout_w7x_standardconfig.nc'
!     n_tolerated_test_arrays_inconsistencies = 0.0
!   
!   geometry_zpinch
!     betaprim = 0.0
!   
!   geometry_from_txt
!     geometry_file = 'input.geometry'
!     overwrite_bmag = .false.
!     overwrite_b_dot_grad_zeta = .false.
!     overwrite_gds2 = .false.
!     overwrite_gds21 = .false.
!     overwrite_gds22 = .false.
!     overwrite_gds23 = .false.
!     overwrite_gds24 = .false.
!     overwrite_gbdrift = .false.
!     overwrite_cvdrift = .false.
!     overwrite_gbdrift0 = .false.
!     set_bmag_const = .false.
! 
! Text options for <geometry_option>:
!    - miller: {default, miller, local}
!    - VMEC: {vmec}
!    - z-pinch: {zpinch}
!    - input profiles: {input.profiles}
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
module namelist_geometry

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_geometry_options
   public :: read_namelist_geometry_miller
   public :: read_namelist_geometry_from_txt
   public :: read_namelist_geometry_vmec
   public :: read_namelist_geometry_zpinch

   ! Parameters need to be public (geometry_option)
   public :: geo_option_local, geo_option_inputprof, geo_option_vmec
   public :: geo_option_multibox, geo_option_zpinch

   ! Parameters need to be public (radial_coordinate)
   public :: radial_coordinate_sgnpsitpsit
   public :: radial_coordinate_minuspsit
   public :: radial_coordinate_r

   private

   ! Create parameters for <geometry_option>
   integer, parameter :: geo_option_local = 1
   integer, parameter :: geo_option_inputprof = 2
   integer, parameter :: geo_option_vmec = 3
   integer, parameter :: geo_option_multibox = 4
   integer, parameter :: geo_option_zpinch = 5

   ! Create parameters for <radial_coordinate>
   integer, parameter :: radial_coordinate_sgnpsitpsit = 1
   integer, parameter :: radial_coordinate_minuspsit = 2
   integer, parameter :: radial_coordinate_r = 3

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                              GEOMETRY OPTIONS                             !
   !****************************************************************************
   subroutine read_namelist_geometry_options(geometry_option_switch, q_as_x)
     
      use mp, only: proc0, mp_abort
      use parameters_physics, only: initialised_parameters_physics
      
      implicit none

      ! Variables that are read from the input file
      integer, intent (out) :: geometry_option_switch
      logical, intent (out) :: q_as_x

      ! Local variable to set <geo_option_switch>
      character (30) :: geometry_option
      
      !-------------------------------------------------------------------------
      
      ! Some geometry flags are based on <radial_variation>
      ! Therefore, we need to read the physics parameters first
      if (.not. initialised_parameters_physics) then
         call mp_abort('Initialise physics parameters before reading geometry namelists. Aborting.')
      end if

      ! Set the default input parameters and read the input file
      if (.not. proc0) return
      call set_default_geometry_options
      call read_input_file_geometry_options
      call check_namelist_geometry_options
      call write_parameters_to_input_file

   contains

      !------------------------ Default input parameters -----------------------
      subroutine set_default_geometry_options
         
         use parameters_physics, only: radial_variation

         implicit none

         ! Text options: {default, miller, local, zpinch, input.profiles, vmec}
         geometry_option = 'local'
         q_as_x = radial_variation

      end subroutine set_default_geometry_options

      !------------------------ Read input file parameters -----------------------
      subroutine read_input_file_geometry_options

         use file_utils, only: input_unit_exist
         use file_units, only: unit_error_file
         use text_options, only: text_option, get_option_value

         implicit none
         
         ! Link text options for <geometry_option> to an integer value
         type(text_option), dimension(6), parameter :: geoopts = (/ &
             text_option('default', geo_option_local), &
             text_option('miller', geo_option_local), &
             text_option('local', geo_option_local), &
             text_option('zpinch', geo_option_zpinch), &
             text_option('input.profiles', geo_option_inputprof), &
             text_option('vmec', geo_option_vmec)/)

         ! Variables in the <geometry_options> namelist
         namelist /geometry_options/ geometry_option, q_as_x
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('geometry_options', dexist)
         if (dexist) read (unit=in_file, nml=geometry_options)

         ! Read the text option in <geometry_option> and store it in <geometry_option_switch>
         call get_option_value(geometry_option, geoopts, geometry_option_switch, &
             unit_error_file, 'geometry_option in namelist_geometry_options.f90')

      end subroutine read_input_file_geometry_options

      !------------------------- Check input parameters ------------------------
      subroutine check_namelist_geometry_options

         use parameters_physics, only: radial_variation
         use file_utils, only: runtype_option_switch, runtype_multibox
         use mp, only: job 

         implicit none
         
         ! Multibox run
         if (radial_variation .and. runtype_option_switch == runtype_multibox .and. job /= 1) then
             geometry_option_switch = geo_option_multibox
         end if

      end subroutine check_namelist_geometry_options
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&geometry_options'
         write (unit, '(A, A, A)') '  geometry_option = "', trim(geometry_option),'"'
         write (unit, '(A, L0)') '  q_as_x = ', q_as_x
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file
     
   end subroutine read_namelist_geometry_options

   !****************************************************************************
   !                    GEOMETRY OPTIONS: GEOMETRY FROM TEXT                   !
   !****************************************************************************

   subroutine read_namelist_geometry_from_txt(geometry_file, overwrite_bmag, &
      overwrite_b_dot_grad_zeta, overwrite_gds2, overwrite_gds21, overwrite_gds22, &
      overwrite_gds23, overwrite_gds24, overwrite_gbdrift, overwrite_cvdrift, &
      overwrite_gbdrift0, set_bmag_const, overwrite_geometry)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      character (100), intent (out) :: geometry_file
      logical, intent (out) :: overwrite_bmag, overwrite_b_dot_grad_zeta, &
         overwrite_gds2, overwrite_gds21, overwrite_gds22, &
         overwrite_gds23, overwrite_gds24, overwrite_gbdrift, &
         overwrite_cvdrift, overwrite_gbdrift0, set_bmag_const
      logical, intent (out) :: overwrite_geometry

      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_geometry_from_txt
      call read_input_file_geometry_from_txt
      call check_namelist_geometry_from_txt 
      call write_parameters_to_input_file

   contains

     !------------------------ Default input parameters -----------------------
     subroutine set_default_geometry_from_txt

         implicit none

         geometry_file = 'input.geometry'
         overwrite_bmag = .false.
         overwrite_b_dot_grad_zeta = .false.
         overwrite_gds2 = .false.
         overwrite_gds21 = .false.
         overwrite_gds22 = .false.
         overwrite_gds23 = .false.
         overwrite_gds24 = .false.
         overwrite_gbdrift = .false.
         overwrite_cvdrift = .false.
         overwrite_gbdrift0 = .false.
         set_bmag_const = .false.
         overwrite_geometry = .false.

     end subroutine set_default_geometry_from_txt

     !------------------------ Read input file parameters -----------------------
     subroutine read_input_file_geometry_from_txt

         use file_utils, only: input_unit_exist

         implicit none

         ! Variables in the <geometry_from_txt> namelist
         namelist /geometry_from_txt/ geometry_file, &
            overwrite_bmag, overwrite_b_dot_grad_zeta, overwrite_gds2, &
            overwrite_gds21, overwrite_gds22, overwrite_gds23, overwrite_gds24, &
            overwrite_cvdrift, overwrite_gbdrift, overwrite_gbdrift0, set_bmag_const

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('geometry_from_txt', dexist)
         if (dexist) read (unit=in_file, nml=geometry_from_txt)

     end subroutine read_input_file_geometry_from_txt

     !------------------------- Check input parameters ------------------------
     subroutine check_namelist_geometry_from_txt

         implicit none

         ! If one of the geometry arrays needs to be overwritten, then <overwrite_geometry> = True
         overwrite_geometry = overwrite_bmag .or. overwrite_b_dot_grad_zeta &
             .or. overwrite_gds2 .or. overwrite_gds21 .or. overwrite_gds22 &
             .or. overwrite_gds23 .or. overwrite_gds24 &
             .or. overwrite_cvdrift .or. overwrite_gbdrift .or. overwrite_gbdrift0

     end subroutine check_namelist_geometry_from_txt
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&geometry_from_txt'
         write (unit, '(A, A, A)') '  geometry_file = "', trim(geometry_file),'"'
         write (unit, '(A, L0)') '  overwrite_bmag = ', overwrite_bmag
         write (unit, '(A, L0)') '  overwrite_b_dot_grad_zeta = ', overwrite_b_dot_grad_zeta
         write (unit, '(A, L0)') '  overwrite_gds2 = ', overwrite_gds2
         write (unit, '(A, L0)') '  overwrite_gds21 = ', overwrite_gds21
         write (unit, '(A, L0)') '  overwrite_gds22 = ', overwrite_gds22
         write (unit, '(A, L0)') '  overwrite_gds23 = ', overwrite_gds23
         write (unit, '(A, L0)') '  overwrite_gds24 = ', overwrite_gds24
         write (unit, '(A, L0)') '  overwrite_cvdrift = ', overwrite_cvdrift
         write (unit, '(A, L0)') '  overwrite_gbdrift = ', overwrite_gbdrift
         write (unit, '(A, L0)') '  overwrite_gbdrift0 = ', overwrite_gbdrift0
         write (unit, '(A, L0)') '  set_bmag_const = ', set_bmag_const
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_geometry_from_txt

   !****************************************************************************
   !                          GEOMETRY OPTIONS: MILLER                         !
   !****************************************************************************
   subroutine read_namelist_geometry_miller(rhoc, rmaj, shift, qinp, shat, &
      kappa, kapprim, tri, triprim, rgeo, betaprim, betadbprim, d2qdr2, d2psidr2, &
      nzed_local, read_profile_variation, write_profile_variation)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      real, intent (out) :: rhoc, rmaj, shift
      real, intent (out) :: qinp, shat, kappa, kapprim
      real, intent (out) :: tri, triprim
      real, intent (out) :: rgeo, betaprim, betadbprim
      real, intent (out) :: d2qdr2, d2psidr2
      integer, intent (out) :: nzed_local
      logical, intent (out) :: write_profile_variation, read_profile_variation

      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_geometry_miller
      call read_input_file_geometry_miller
      call write_parameters_to_input_file
      
   contains

     !------------------------ Default input parameters -----------------------
     subroutine set_default_geometry_miller

         implicit none

         rhoc = 0.5
         rmaj = 2.77778
         shift = 0.0
         qinp = 1.4
         shat = 0.796
         kappa = 1.0
         kapprim = 0.0
         tri = 0.0
         triprim = 0.0
         rgeo = 2.77778
         betaprim = 0.0       ! betaprim = -(4pi/Bref^2)*d(ptot)/drho
         betadbprim = 0.0     ! betadbprim = -(4pi/Bref^2)*d^2ptot/drho^2
         d2qdr2 = 0.0
         d2psidr2 = 0.0
         nzed_local = 128
         read_profile_variation = .false.
         write_profile_variation = .false.

     end subroutine set_default_geometry_miller

     !------------------------ Read input file parameters -----------------------
     subroutine read_input_file_geometry_miller

         use file_utils, only: input_unit_exist

         implicit none

         ! Variables in the <geometry_miller> namelist
         namelist /geometry_miller/ rhoc, rmaj, shift, qinp, shat, &
            kappa, kapprim, tri, triprim, rgeo, betaprim, betadbprim, d2qdr2, &
            d2psidr2, nzed_local, read_profile_variation, write_profile_variation

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('geometry_miller', dexist)
         if (dexist) read (unit=in_file, nml=geometry_miller)

     end subroutine read_input_file_geometry_miller
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&geometry_miller'
         write (unit, '(A, ES0.4)') '  rhoc = ', rhoc
         write (unit, '(A, ES0.4)') '  rmaj = ', rmaj
         write (unit, '(A, ES0.4)') '  shift = ', shift
         write (unit, '(A, ES0.4)') '  qinp = ', qinp
         write (unit, '(A, ES0.4)') '  shat = ', shat
         write (unit, '(A, ES0.4)') '  kappa = ', kappa
         write (unit, '(A, ES0.4)') '  kapprim = ', kapprim
         write (unit, '(A, ES0.4)') '  tri = ', tri
         write (unit, '(A, ES0.4)') '  triprim = ', triprim
         write (unit, '(A, ES0.4)') '  rgeo = ', rgeo
         write (unit, '(A, ES0.4)') '  betaprim = ', betaprim
         write (unit, '(A, ES0.4)') '  betadbprim = ', betadbprim
         write (unit, '(A, ES0.4)') '  d2qdr2 = ', d2qdr2
         write (unit, '(A, ES0.4)') '  d2psidr2 = ', d2psidr2
         write (unit, '(A, I0)') '  nzed_local = ', nzed_local
         write (unit, '(A, L0)') '  read_profile_variation = ', read_profile_variation
         write (unit, '(A, L0)') '  write_profile_variation = ', write_profile_variation
         write (unit, '(A)') '/'
         write (unit, '(A)') ''
      
      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_geometry_miller

   !****************************************************************************
   !                          GEOMETRY OPTIONS: VMEC                           !
   !****************************************************************************
   subroutine read_namelist_geometry_vmec(alpha0, zeta_center, rectangular_cross_section, & 
      nfield_periods, torflux, z_grid_refinement_factor, surface_option, radial_coordinate_switch, &
      verbose, vmec_filename, n_tolerated_test_arrays_inconsistencies)

      use mp, only: proc0, mp_abort
      use grids_z, only: initialised_grids_z

      implicit none

      ! Variables that are read from the input file
      real, intent (out) :: alpha0, zeta_center
      real, intent (out) :: nfield_periods, torflux
      logical, intent (out) :: verbose, rectangular_cross_section
      integer, intent (out) :: radial_coordinate_switch
      integer, intent (out) :: n_tolerated_test_arrays_inconsistencies
      integer, intent (out) :: z_grid_refinement_factor
      integer, intent (out) :: surface_option
      character(2000), intent (out) :: vmec_filename

      ! Local variable to set <radial_coordinate_switch>
      character(20) :: radial_coordinate
      
      !-------------------------------------------------------------------------
      
      ! The geometry options need <zed_equal_arc> from the grids_z module
      ! Therefore, we need to read the physics parameters first
      if (.not. initialised_grids_z) then
         call mp_abort('Initialise z-grid parameters before reading geometry namelists. Aborting.')
      end if

      ! Set the default input parameters and read the input file
      if (.not. proc0) return
      call set_default_geometry_vmec
      call read_input_file_geometry_vmec
      call write_parameters_to_input_file

   contains
   
      !------------------------ Default input parameters -----------------------
      subroutine set_default_geometry_vmec

         use grids_z, only: zed_equal_arc

         implicit none

         ! Default parameters for VMEC
         vmec_filename = 'equilibria/wout_w7x_standardConfig.nc'
         alpha0 = 0.0
         zeta_center = 0.0
         nfield_periods = -1.0
         torflux = 0.6354167d+0
         surface_option = 0
         verbose = .true.
         n_tolerated_test_arrays_inconsistencies = 0
         z_grid_refinement_factor = 1
         radial_coordinate = 'sgn(psi_t) psi_t'

         ! If we use the normalized arc-length as the parallel coordinate,
         ! then use <z_grid_refinement_factor> more z-points to calculate
         ! the geometry arrays in VMEC, by using a smaller step dzeta.
         if (zed_equal_arc) then
             z_grid_refinement_factor = 4
         else
             z_grid_refinement_factor = 1
         end if

         ! For alpha=0 the perpendicular cross-section is rectangular at zeta_center=0
         ! For alpha!=0 it is not in the original stella, since alpha = theta_p - iota*zeta
         ! To make the cross-section rectangular we need to define alpa = theta_p - iota(zeta-zeta_center) 
         ! This is done by toggling <rectangular_cross_section> to .true.
         rectangular_cross_section = .false.

      end subroutine set_default_geometry_vmec

     !------------------------ Read input file parameters -----------------------
     subroutine read_input_file_geometry_vmec

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr
         
         ! Link text options for <radial_coordinate> to an integer value
         type(text_option), dimension(5), parameter :: radial_coordinate_options = & 
             (/text_option('default', radial_coordinate_sgnpsitpsit), &
               text_option('sgn(psi_t) psi_t', radial_coordinate_sgnpsitpsit), &
               text_option('original', radial_coordinate_minuspsit), &
               text_option('- psi_t', radial_coordinate_minuspsit), &
               text_option('r', radial_coordinate_r)/)

         ! Variables in the <geometry_vmec> namelist
         namelist /geometry_vmec/ alpha0, zeta_center, rectangular_cross_section, nfield_periods, &
             torflux, z_grid_refinement_factor, surface_option, radial_coordinate, &
             verbose, vmec_filename, n_tolerated_test_arrays_inconsistencies
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('geometry_vmec', dexist)
         if (dexist) read (unit=in_file, nml=geometry_vmec)

         ! Read the text option in <radial_coordinate> and store it in <radial_coordinate_switch>
         ierr = error_unit()
         call get_option_value(radial_coordinate, radial_coordinate_options, &
            radial_coordinate_switch, ierr, 'radial_coordinate in namelist_geometry.f90')

     end subroutine read_input_file_geometry_vmec
     
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&geometry_vmec'
         write (unit, '(A, A, A)') '  vmec_filename = "', trim(vmec_filename),'"'
         write (unit, '(A, ES0.4)') '  alpha0 = ', alpha0
         write (unit, '(A, ES0.4)') '  zeta_center = ', zeta_center
         write (unit, '(A, ES0.8)') '  nfield_periods = ', nfield_periods
         write (unit, '(A, ES0.4)') '  torflux = ', torflux
         write (unit, '(A, I0)') '  z_grid_refinement_factor = ', z_grid_refinement_factor 
         write (unit, '(A, I0)') '  surface_option = ', surface_option 
         write (unit, '(A, I0)') '  n_tolerated_test_arrays_inconsistencies = ', n_tolerated_test_arrays_inconsistencies
         write (unit, '(A, L0)') '  verbose = ', verbose
         write (unit, '(A, L0)') '  rectangular_cross_section = ', rectangular_cross_section
         write (unit, '(A, A, A)') '  radial_coordinate = "', trim(radial_coordinate),'"'
         write (unit, '(A)') '/'
         write (unit, '(A)') ''
      
      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_geometry_vmec

   !****************************************************************************
   !                         GEOMETRY OPTIONS: Z-PINCH                         !
   !****************************************************************************
   subroutine read_namelist_geometry_zpinch(betaprim)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      real, intent (out) :: betaprim

      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_geometry_zpinch
      call read_input_file_geometry_zpinch
      call write_parameters_to_input_file

   contains

      !------------------------ Default input parameters -----------------------
      subroutine set_default_geometry_zpinch

         implicit none

         betaprim = 0.0

      end subroutine set_default_geometry_zpinch

      !------------------------ Read input file parameters -----------------------
      subroutine read_input_file_geometry_zpinch

         use file_utils, only: input_unit_exist

         implicit none

         namelist /geometry_zpinch/ betaprim
         in_file = input_unit_exist('geometry_zpinch', dexist)
         if (dexist) read (unit=in_file, nml=geometry_zpinch)

      end subroutine read_input_file_geometry_zpinch
     
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&geometry_zpinch'
         write (unit, '(A, ES0.4)') '  betaprim = ', betaprim
         write (unit, '(A)') '/'
         write (unit, '(A)') ''
      
      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_geometry_zpinch

end module namelist_geometry
