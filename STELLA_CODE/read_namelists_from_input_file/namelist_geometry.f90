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
!     rmaj = 2.77778
!     shift = 0.0
!     qinp = 1.4
!     shat = 0.796
!     kappa = 1.0
!     kapprim = 0.0
!     tri = 0.0
!     triprim = 0.0
!     rgeo = 2.77778
!     betaprim = 0.0
!     betadbprim = 0.0
!     d2qdr2 = 0.0
!     d2psidr2 = 0.0
!     nzed_local = 128
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
!     surface_option = 0
!     radial_coordinate = 'sgn(psi_t)psi_t'
!     verbose = .true.
!     vmec_filename = 'wout_w7x_standardConfig.nc'
!     n_tolerated_test_arrays_inconsistencies = 0
!   
!   geometry_zpinch
!     betaprim = 0.0
!   
!   geometry_from_txt
!     geometry_file = 'input.geometry'
!     overwrite_bmag = .false.
!     overwrite_b_dot_gradzeta = .false.
!     overwrite_grady_dot_grady = .false.
!     overwrite_gradx_dot_grady = .false.
!     overwrite_gradx_dot_gradx = .false.
!     overwrite_gds23 = .false.
!     overwrite_gds24 = .false.
!     overwrite_B_times_gradB_dot_grady = .false.
!     overwrite_B_times_kappa_dot_grady = .false.
!     overwrite_B_times_gradB_dot_gradx = .false.
!
!   <geometry_file> must be a *.geometry file previously written by stella
!   itself (write_geometric_coefficients in geometry.f90); overwrite_gds23
!   and overwrite_gds24 abort at runtime because gds23/gds24 are no longer
!   among the columns that file contains.
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
   public :: read_namelist_geometry_gvec
   public :: read_namelist_geometry_zpinch
   public :: read_namelist_geometry_slab

   ! Parameters need to be public (geometry_option)
   public :: geo_option_local, geo_option_inputprof, geo_option_vmec
   public :: geo_option_multibox, geo_option_zpinch, geo_option_slab, geo_option_gvec

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
   integer, parameter :: geo_option_slab = 6
   integer, parameter :: geo_option_gvec = 7

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
         type(text_option), dimension(8), parameter :: geoopts = (/ &
             text_option('default', geo_option_local), &
             text_option('miller', geo_option_local), &
             text_option('local', geo_option_local), &
             text_option('zpinch', geo_option_zpinch), &
             text_option('input.profiles', geo_option_inputprof), &
             text_option('vmec', geo_option_vmec), &
             text_option('slab', geo_option_slab), &
             text_option('gvec', geo_option_gvec)/)

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
         use grids_z, only: nperiod
         use mp, only: job, mp_abort

         implicit none
         
         ! Abort if nperiod > 1 for VMEC or GVEC geometries
         if (geometry_option_switch==geo_option_vmec .or. &
             geometry_option_switch==geo_option_gvec) then
            if (nperiod/=1) then
               call mp_abort('If geometry_option = "vmec" or "gvec" then nperiod should be set to 1. Aborting.')
            end if
         end if
            
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
         write (unit, '(A, L1)') '  q_as_x = ', q_as_x
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file
     
   end subroutine read_namelist_geometry_options

   !****************************************************************************
   !                    GEOMETRY OPTIONS: GEOMETRY FROM TEXT                   !
   !****************************************************************************

   subroutine read_namelist_geometry_from_txt(geometry_file, overwrite_bmag, &
      overwrite_b_dot_gradzeta, overwrite_grady_dot_grady, overwrite_gradx_dot_grady, overwrite_gradx_dot_gradx, &
      overwrite_gds23, overwrite_gds24, overwrite_B_times_gradB_dot_grady, overwrite_B_times_kappa_dot_grady, &
      overwrite_B_times_gradB_dot_gradx, overwrite_geometry)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      character (100), intent (out) :: geometry_file
      logical, intent (out) :: overwrite_bmag, overwrite_b_dot_gradzeta, &
         overwrite_grady_dot_grady, overwrite_gradx_dot_grady, overwrite_gradx_dot_gradx, &
         overwrite_gds23, overwrite_gds24, overwrite_B_times_gradB_dot_grady, &
         overwrite_B_times_kappa_dot_grady, overwrite_B_times_gradB_dot_gradx
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
         overwrite_b_dot_gradzeta = .false.
         overwrite_grady_dot_grady = .false.
         overwrite_gradx_dot_grady = .false.
         overwrite_gradx_dot_gradx = .false.
         overwrite_gds23 = .false.
         overwrite_gds24 = .false.
         overwrite_B_times_gradB_dot_grady = .false.
         overwrite_B_times_kappa_dot_grady = .false.
         overwrite_B_times_gradB_dot_gradx = .false.
         overwrite_geometry = .false.

     end subroutine set_default_geometry_from_txt

     !------------------------ Read input file parameters -----------------------
     subroutine read_input_file_geometry_from_txt

         use file_utils, only: input_unit_exist

         implicit none

         ! Variables in the <geometry_from_txt> namelist
         namelist /geometry_from_txt/ geometry_file, &
            overwrite_bmag, overwrite_b_dot_gradzeta, overwrite_grady_dot_grady, &
            overwrite_gradx_dot_grady, overwrite_gradx_dot_gradx, overwrite_gds23, overwrite_gds24, &
            overwrite_B_times_kappa_dot_grady, overwrite_B_times_gradB_dot_grady, overwrite_B_times_gradB_dot_gradx

         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist('geometry_from_txt', dexist)
         if (dexist) read (unit=in_file, nml=geometry_from_txt)

     end subroutine read_input_file_geometry_from_txt

     !------------------------- Check input parameters ------------------------
     subroutine check_namelist_geometry_from_txt

         implicit none

         ! If one of the geometry arrays needs to be overwritten, then <overwrite_geometry> = True
         overwrite_geometry = overwrite_bmag .or. overwrite_b_dot_gradzeta &
             .or. overwrite_grady_dot_grady .or. overwrite_gradx_dot_grady .or. overwrite_gradx_dot_gradx &
             .or. overwrite_gds23 .or. overwrite_gds24 &
             .or. overwrite_B_times_kappa_dot_grady .or. overwrite_B_times_gradB_dot_grady .or. overwrite_B_times_gradB_dot_gradx

     end subroutine check_namelist_geometry_from_txt
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&geometry_from_txt'
         write (unit, '(A, A, A)') '  geometry_file = "', trim(geometry_file),'"'
         write (unit, '(A, L1)') '  overwrite_bmag = ', overwrite_bmag
         write (unit, '(A, L1)') '  overwrite_b_dot_gradzeta = ', overwrite_b_dot_gradzeta
         write (unit, '(A, L1)') '  overwrite_grady_dot_grady = ', overwrite_grady_dot_grady
         write (unit, '(A, L1)') '  overwrite_gradx_dot_grady = ', overwrite_gradx_dot_grady
         write (unit, '(A, L1)') '  overwrite_gradx_dot_gradx = ', overwrite_gradx_dot_gradx
         write (unit, '(A, L1)') '  overwrite_gds23 = ', overwrite_gds23
         write (unit, '(A, L1)') '  overwrite_gds24 = ', overwrite_gds24
         write (unit, '(A, L1)') '  overwrite_B_times_kappa_dot_grady = ', overwrite_B_times_kappa_dot_grady
         write (unit, '(A, L1)') '  overwrite_B_times_gradB_dot_grady = ', overwrite_B_times_gradB_dot_grady
         write (unit, '(A, L1)') '  overwrite_B_times_gradB_dot_gradx = ', overwrite_B_times_gradB_dot_gradx
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
         write (unit, '(A, ES12.4)') '  rhoc = ', rhoc
         write (unit, '(A, ES12.4)') '  rmaj = ', rmaj
         write (unit, '(A, ES12.4)') '  shift = ', shift
         write (unit, '(A, ES12.4)') '  qinp = ', qinp
         write (unit, '(A, ES12.4)') '  shat = ', shat
         write (unit, '(A, ES12.4)') '  kappa = ', kappa
         write (unit, '(A, ES12.4)') '  kapprim = ', kapprim
         write (unit, '(A, ES12.4)') '  tri = ', tri
         write (unit, '(A, ES12.4)') '  triprim = ', triprim
         write (unit, '(A, ES12.4)') '  rgeo = ', rgeo
         write (unit, '(A, ES12.4)') '  betaprim = ', betaprim
         write (unit, '(A, ES12.4)') '  betadbprim = ', betadbprim
         write (unit, '(A, ES12.4)') '  d2qdr2 = ', d2qdr2
         write (unit, '(A, ES12.4)') '  d2psidr2 = ', d2psidr2
         write (unit, '(A, I0)') '  nzed_local = ', nzed_local
         write (unit, '(A, L1)') '  read_profile_variation = ', read_profile_variation
         write (unit, '(A, L1)') '  write_profile_variation = ', write_profile_variation
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
         vmec_filename = 'wout_w7x_standardConfig.nc'
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
         write (unit, '(A, ES12.4)') '  alpha0 = ', alpha0
         write (unit, '(A, ES12.4)') '  zeta_center = ', zeta_center
         write (unit, '(A, ES16.8)') '  nfield_periods = ', nfield_periods
         write (unit, '(A, ES12.4)') '  torflux = ', torflux
         write (unit, '(A, I0)') '  z_grid_refinement_factor = ', z_grid_refinement_factor 
         write (unit, '(A, I0)') '  surface_option = ', surface_option 
         write (unit, '(A, I0)') '  n_tolerated_test_arrays_inconsistencies = ', n_tolerated_test_arrays_inconsistencies
         write (unit, '(A, L1)') '  verbose = ', verbose
         write (unit, '(A, L1)') '  rectangular_cross_section = ', rectangular_cross_section
         write (unit, '(A, A, A)') '  radial_coordinate = "', trim(radial_coordinate),'"'
         write (unit, '(A)') '/'
         write (unit, '(A)') ''
      
      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_geometry_vmec

   !****************************************************************************
   !                           GEOMETRY OPTIONS: GVEC                         !
   !****************************************************************************
   ! Reads the &geometry_gvec namelist.
   !
   ! Parameters mirror the &geometry_vmec namelist:
   !   gvec_filename            - netCDF file from gvec_to_stella.py
   !   torflux                  - normalized toroidal flux s ∈ (0,1] of the surface
   !   alpha0                   - label of the first field line (radians)
   !   zeta_center              - centre of the toroidal domain (radians)
   !   nfield_periods           - number of field periods to follow (<= nfp; -1 = all)
   !   z_grid_refinement_factor - refinement relative to stella's nzgrid
   !                              (file must have been generated with
   !                               --nzgrid_file = nzgrid * z_grid_refinement_factor)
   !   verbose                  - print geometry information
   !
   ! The grid parameters (torflux, alpha0, zeta_center, nfield_periods) are
   ! validated against the pre-computed netCDF file.  If they do not match, stella
   ! aborts with an error message showing which gvec_to_stella.py arguments to use.
   !****************************************************************************
   subroutine read_namelist_geometry_gvec( &
      z_grid_refinement_factor, radial_coordinate_switch, verbose, gvec_filename, &
      torflux, alpha0, zeta_center, nfield_periods, &
      gvec_state_file, gvec_parameter_file, sfl_coord, sfl_factor)

      use mp,      only: proc0, mp_abort
      use grids_z, only: initialised_grids_z, zed_equal_arc

      implicit none

      integer,         intent(out) :: z_grid_refinement_factor
      integer,         intent(out) :: radial_coordinate_switch
      logical,         intent(out) :: verbose
      character(2000), intent(out) :: gvec_filename
      real,            intent(out) :: torflux, alpha0, zeta_center, nfield_periods
      ! Direct path: read a GVEC state file through the GVEC library instead of
      ! a geometry netCDF written by gvec_to_stella.py.
      character(2000), intent(out) :: gvec_state_file, gvec_parameter_file
      integer,         intent(out) :: sfl_coord, sfl_factor

      !-------------------------------------------------------------------------

      if (.not. initialised_grids_z) then
         call mp_abort('Initialise z-grid parameters before reading geometry namelists. Aborting.')
      end if

      if (.not. proc0) return
      call set_defaults
      call read_from_file
      call write_to_output

   contains

      subroutine set_defaults
         implicit none
         gvec_filename            = 'gvec_stella_geometry.nc'
         gvec_state_file          = ''      ! empty -> use the netCDF path
         gvec_parameter_file      = ''      ! needed for non-cylindrical h-maps
         sfl_coord                = 1       ! 0: old PEST, 1: PEST, 2: Boozer
         sfl_factor               = 4       ! Fourier-mode factor for the SFL transform
         verbose                  = .true.
         radial_coordinate_switch = radial_coordinate_sgnpsitpsit
         torflux                  = -1.0   ! -1 → read from file
         alpha0                   = -999.0 ! -999 → read from file
         zeta_center              = -999.0 ! -999 → read from file
         nfield_periods           = -1.0   ! -1   → read from file (all periods)
         if (zed_equal_arc) then
            z_grid_refinement_factor = 4
         else
            z_grid_refinement_factor = 1
         end if
      end subroutine set_defaults

      subroutine read_from_file
         use file_utils, only: input_unit_exist
         implicit none
         namelist /geometry_gvec/ gvec_filename, verbose, z_grid_refinement_factor, &
            torflux, alpha0, zeta_center, nfield_periods, &
            gvec_state_file, gvec_parameter_file, sfl_coord, sfl_factor
         in_file = input_unit_exist('geometry_gvec', dexist)
         if (dexist) read (unit=in_file, nml=geometry_gvec)
      end subroutine read_from_file

      subroutine write_to_output
         use file_units, only: unit => unit_input_file_with_defaults
         implicit none
         write (unit, '(A)')        '&geometry_gvec'
         write (unit, '(2A)')       '  gvec_filename = "', trim(gvec_filename)//'"'
         write (unit, '(2A)')       '  gvec_state_file = "', trim(gvec_state_file)//'"'
         write (unit, '(2A)')       '  gvec_parameter_file = "', trim(gvec_parameter_file)//'"'
         write (unit, '(A,I0)')     '  sfl_coord = ', sfl_coord
         write (unit, '(A,I0)')     '  sfl_factor = ', sfl_factor
         write (unit, '(A,ES10.4)') '  torflux = ', torflux
         write (unit, '(A,ES10.4)') '  alpha0 = ', alpha0
         write (unit, '(A,ES10.4)') '  zeta_center = ', zeta_center
         write (unit, '(A,ES10.4)') '  nfield_periods = ', nfield_periods
         write (unit, '(A,I0)')     '  z_grid_refinement_factor = ', z_grid_refinement_factor
         write (unit, '(A,L1)')     '  verbose = ', verbose
         write (unit, '(A)')        '/'
         write (unit, '(A)')        ''
      end subroutine write_to_output

   end subroutine read_namelist_geometry_gvec

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
         write (unit, '(A, ES12.4)') '  betaprim = ', betaprim
         write (unit, '(A)') '/'
         write (unit, '(A)') ''
      
      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_geometry_zpinch

   !****************************************************************************
   !                          GEOMETRY OPTIONS: SLAB                           !
   !****************************************************************************
   subroutine read_namelist_geometry_slab(shat, qinp, rmaj, betaprim, curvy_slab)

      use mp, only: proc0, mp_abort

      implicit none

      ! Variables that are read from the input file
      real, intent(out) :: shat, qinp, rmaj, betaprim
      logical, intent(out) :: curvy_slab

      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_geometry_slab
      call read_input_file_geometry_slab
      call check_namelist_geometry_slab
      call write_parameters_to_input_file

   contains

      !------------------------ Default input parameters -----------------------
      subroutine set_default_geometry_slab

         implicit none

         shat = 0.0
         qinp = 1.4
         rmaj = 3.0
         betaprim = 0.0
         curvy_slab = .false.

      end subroutine set_default_geometry_slab

      !------------------------ Read input file parameters -----------------------
      subroutine read_input_file_geometry_slab

         use file_utils, only: input_unit_exist

         implicit none

         namelist /geometry_slab/ shat, qinp, rmaj, betaprim, curvy_slab
         in_file = input_unit_exist('geometry_slab', dexist)
         if (dexist) read (unit=in_file, nml=geometry_slab)

      end subroutine read_input_file_geometry_slab

      !------------------------- Check input parameters ------------------------
      subroutine check_namelist_geometry_slab

         implicit none

         if (qinp <= 0.0) then
            call mp_abort('geometry_slab: qinp must be positive. Aborting.')
         end if
         if (rmaj <= 0.0) then
            call mp_abort('geometry_slab: rmaj must be positive. Aborting.')
         end if

      end subroutine check_namelist_geometry_slab

      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         write (unit, '(A)') '&geometry_slab'
         write (unit, '(A, ES12.4)') '  shat = ', shat
         write (unit, '(A, ES12.4)') '  qinp = ', qinp
         write (unit, '(A, ES12.4)') '  rmaj = ', rmaj
         write (unit, '(A, ES12.4)') '  betaprim = ', betaprim
         write (unit, '(A, L1)') '  curvy_slab = ', curvy_slab
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_geometry_slab

end module namelist_geometry
