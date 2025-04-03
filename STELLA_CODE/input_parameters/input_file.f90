!###############################################################################
!########################## READ ALL STELLA NAMELISTS ##########################
!###############################################################################
! 
! This module will set the default input input parameters for each name list,
! and it will read the stella input file per namelist.
! 
! For each namelists two routines will exist:
!    - read_namelist_<namelist>
!    - set_maxwellian_parameters_<namelist>
! 
! First we will set the default input parameters, and then we will overwrite 
! any default options with those specified in the input file. Optionally
! we can check if any input variables are clashing with each other.
! 
! Overview of stella namelists:
! 
! GEOMERTRY
!   geometry_option
!   overwrite_geometry
!   geometry_vmec (renamed from vmec_parameters)
!   geometry_miller (renamed from millergeo_parameters)
! 
! PHYSICS
!   gyrokinetic_terms
!   scale_gyrokinetic_terms
!   adiabatic_electron_response
!   adiabatic_ion_response
!   electromagnetic
!   full_flux_surface
!   extra_physics 
! 
! GRIDS
!   vpamu_grid
!   z_grid (renamed from zgrid_parameters)
!   z_boundary_condition
!   species_knobs
!   species_parameters_1
!   species_parameters_2
!   kxky_grid_option
!   kxky_grid_range
!   kxky_grid_box
! 
! DIAGNOSTICS
!   diagnostics
!   diagnostics_potential
!   diagnostics_omega
!   diagnostics_distribution
!   diagnostics_fluxes
!   diagnostics_moments
! 
! INITIALIZE FIELDS
!   initialize_distribution (renamed from init_g_knobs)
!   initialize_distribution_maxwellian
!   initialize_distribution_noise
!   initialize_distribution_kpar
!   initialize_distribution_rh
!   restart_options
! 
! DISSIPATION AND COLLISIONS
!   dissipation_and_collisions_options (Renamed from &dissipation)
!   collisions_dougherty
!   collisions_fokker_planck (Renamed from &collisions_fp)
!   hyper_dissipation
! 
! TIME TRACE
!   time_trace_options
!   time_step
! 
! NUMERICS
!   numerical_algorithms
!   numerical_upwinding_for_derivatives
! 
! NEOCLASSICS
!   neoclassical_input
!   euterpe_parameters
!   sources
! 
! RADIAL VARIATION
!   radial_variation
! 
! PARALLELISATION
!   parallelisation (renamed from &layout_knobs)
! 
! VERBOSE
!   debug_flags
! 
!###############################################################################
module input_file

   implicit none

   public :: read_namelist_dissipation
   public :: read_namelist_initialize_distribution
   public :: read_namelist_initialize_distribution_maxwellian
   public :: read_namelist_initialize_distribution_noise
   public :: read_namelist_initialize_distribution_kpar
   public :: read_namelist_initialize_distribution_rh
   public :: read_namelist_restart_options
   
   ! Parameters need to be public
   public :: init_distribution_option_maxwellian, init_distribution_option_noise, init_distribution_option_restart_many
   public :: init_distribution_option_kpar, init_distribution_option_rh, init_distribution_option_remap
   
   private
   
   ! Create parameters for the <init_distribution_option>
   integer, parameter :: init_distribution_option_maxwellian = 1
   integer, parameter :: init_distribution_option_noise = 2
   integer, parameter :: init_distribution_option_restart_many = 3
   integer, parameter :: init_distribution_option_kpar = 4
   integer, parameter :: init_distribution_option_rh = 5
   integer, parameter :: init_distribution_option_remap = 6 
   
   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   !****************************************************************************
   !                                DISSIPATION                                !
   !****************************************************************************
   subroutine read_namelist_dissipation(include_collisions, collisions_implicit, collision_model, hyper_dissipation)

      use mp, only: proc0

      implicit none

      logical, intent(out) :: include_collisions, collisions_implicit, hyper_dissipation
      character(30), intent(out) :: collision_model

      if (.not. proc0) return
      call set_maxwellian_parameters_dissipation
      call read_input_file_dissipation
      call check_inputs_dissipation

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_maxwellian_parameters_dissipation

         implicit none

         ! By default we do not include collisions nor dissipation
         include_collisions = .false.
         collisions_implicit = .true.
         hyper_dissipation = .false.

         ! Options: dougherty or fokker-planck
         collision_model = "dougherty"

      end subroutine set_maxwellian_parameters_dissipation

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_dissipation

         use file_utils, only: input_unit_exist
         implicit none

         namelist /dissipation/ include_collisions, collisions_implicit, collision_model, hyper_dissipation
         in_file = input_unit_exist("dissipation", dexist)
         if (dexist) read (unit=in_file, nml=dissipation)

      end subroutine read_input_file_dissipation

      !------------------------- Check input parameters ------------------------
      subroutine check_inputs_dissipation

         implicit none

         if (.not. include_collisions) collisions_implicit = .false.

      end subroutine check_inputs_dissipation

   end subroutine read_namelist_dissipation
   
   !****************************************************************************
   !                   INITIALIZE POTENTIAL: READ THE SWITCH                   !
   !****************************************************************************
   subroutine read_namelist_initialize_distribution(init_distribution_switch, phiinit, scale_to_phiinit)

      use mp, only: proc0

      implicit none
      
      ! Variables that are read from the input file
      real, intent(out) :: phiinit
      logical, intent(out) :: scale_to_phiinit
      integer, intent(out) :: init_distribution_switch
      
      ! Local variable to set <init_distribution_switch>
      character(20) :: initialize_distribution_option

      if (.not. proc0) return
      call set_maxwellian_parameters_initialize_distribution
      call read_input_file_initialize_distribution

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_maxwellian_parameters_initialize_distribution

         implicit none

         ! Options: {default, maxwellian, snoise, many, kpar, rh, remap}
         initialize_distribution_option = "maxwellian"
         
         ! Other options
         phiinit = 1.0
         scale_to_phiinit = .false.

      end subroutine set_maxwellian_parameters_initialize_distribution

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_initialize_distribution

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr
      
         ! Link text options for <initialize_distribution_option> to an integer value
         type(text_option), dimension(7), parameter :: init_distribution_options = &
             (/text_option('default', init_distribution_option_maxwellian), &
               text_option('maxwellian', init_distribution_option_maxwellian), &
               text_option('noise', init_distribution_option_noise), &
               text_option('many', init_distribution_option_restart_many), &
               text_option('kpar', init_distribution_option_kpar), &
               text_option('rh', init_distribution_option_rh), &
               text_option('remap', init_distribution_option_remap)/)

         ! Variables in the <initialize_distribution> namelist
         namelist /initialize_distribution/ initialize_distribution_option, phiinit, scale_to_phiinit
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist("initialize_distribution", dexist)
         if (dexist) read (unit=in_file, nml=initialize_distribution)
         
         ! Read the text option in <initialize_distribution> and store it in <init_distribution_switch>
         ierr = error_unit()
         call get_option_value(initialize_distribution_option, init_distribution_options, init_distribution_switch, &
            ierr, "initialize_distribution_option in initialize_distribution")

      end subroutine read_input_file_initialize_distribution

   end subroutine read_namelist_initialize_distribution
   
   !****************************************************************************
   !                      INITIALIZE POTENTIAL: MAXWELLIAN                     !
   !****************************************************************************
   subroutine read_namelist_initialize_distribution_maxwellian(width0, den0, upar0, oddparity, left, chop_side)

      use mp, only: proc0

      implicit none
      
      real, intent(out) :: width0, den0, upar0
      logical, intent(out) :: oddparity
      logical, intent(out) :: left
      logical, intent(out) :: chop_side

      if (.not. proc0) return
      call set_maxwellian_parameters_initialize_distribution_maxwellian
      call read_input_file_initialize_distribution_maxwellian

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_maxwellian_parameters_initialize_distribution_maxwellian

         implicit none
         
         width0 = -3.5
         den0 = 1.
         upar0 = 0.
         oddparity = .false.
         chop_side = .false.
         left = .true.

      end subroutine set_maxwellian_parameters_initialize_distribution_maxwellian

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_initialize_distribution_maxwellian

         use file_utils, only: input_unit_exist
         implicit none

         namelist /initialize_distribution_maxwellian/ width0, den0, upar0, oddparity, left, chop_side
         in_file = input_unit_exist("initialize_distribution_maxwellian", dexist)
         if (dexist) read (unit=in_file, nml=initialize_distribution_maxwellian) 

      end subroutine read_input_file_initialize_distribution_maxwellian

   end subroutine read_namelist_initialize_distribution_maxwellian
   
   !****************************************************************************
   !                       INITIALIZE DISTRIBUTION: NOISE                      !
   !****************************************************************************
   subroutine read_namelist_initialize_distribution_noise(zf_init, left, chop_side)

      use mp, only: proc0

      implicit none
      
      real, intent(out) :: zf_init
      logical, intent(out) :: left
      logical, intent(out) :: chop_side

      if (.not. proc0) return
      call set_default_parameters_initialize_distribution_noise
      call read_input_file_initialize_distribution_noise

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_initialize_distribution_noise

         implicit none
         
         zf_init = 1.0
         chop_side = .false.
         left = .true.

      end subroutine set_default_parameters_initialize_distribution_noise

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_initialize_distribution_noise

         use file_utils, only: input_unit_exist
         implicit none

         namelist /initialize_distribution_noise/ zf_init, left, chop_side
         in_file = input_unit_exist("initialize_distribution_noise", dexist)
         if (dexist) read (unit=in_file, nml=initialize_distribution_noise) 

      end subroutine read_input_file_initialize_distribution_noise

   end subroutine read_namelist_initialize_distribution_noise
   
   !****************************************************************************
   !                         INITIALIZE POTENTIAL: KPAR                        !
   !****************************************************************************
   subroutine read_namelist_initialize_distribution_kpar(&
         width0, refac, imfac, den0, upar0, tpar0, tperp0, &
         den1, upar1, tpar1, tperp1, den2, upar2, tpar2, tperp2, left, chop_side)

      use mp, only: proc0

      implicit none
      
      real, intent(out) :: width0, imfac, refac
      real, intent(out) :: den0, upar0, tpar0, tperp0
      real, intent(out) :: den1, upar1, tpar1, tperp1
      real, intent(out) :: den2, upar2, tpar2, tperp2
      logical, intent(out) :: left
      logical, intent(out) :: chop_side

      if (.not. proc0) return
      call set_default_parameters_initialize_distribution_kpar
      call read_input_file_initialize_distribution_kpar

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_initialize_distribution_kpar

         implicit none
         
         width0 = -3.5
         refac = 1.
         imfac = 0.
         den0 = 1.
         upar0 = 0.
         tpar0 = 0.
         tperp0 = 0.
         den1 = 0.
         upar1 = 0.
         tpar1 = 0.
         tperp1 = 0.
         den2 = 0.
         upar2 = 0.
         tpar2 = 0.
         tperp2 = 0.
         chop_side = .false.
         left = .true.

      end subroutine set_default_parameters_initialize_distribution_kpar

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_initialize_distribution_kpar

         use file_utils, only: input_unit_exist
         implicit none

         namelist /initialize_distribution_kpar/ width0, refac, imfac, den0, upar0, &
            tpar0, tperp0, den1, upar1, tpar1, tperp1, den2, upar2, tpar2, tperp2, left, chop_side
         in_file = input_unit_exist("initialize_distribution_kpar", dexist)
         if (dexist) read (unit=in_file, nml=initialize_distribution_kpar) 

      end subroutine read_input_file_initialize_distribution_kpar

   end subroutine read_namelist_initialize_distribution_kpar
   
   !****************************************************************************
   !                          INITIALIZE POTENTIAL: RH                         !
   !****************************************************************************
   subroutine read_namelist_initialize_distribution_rh(kxmin, kxmax, imfac, refac)

      use mp, only: proc0

      implicit none
      
      real, intent(out) :: kxmax, kxmin
      real, intent(out) :: imfac, refac

      if (.not. proc0) return
      call set_default_parameters_initialize_distribution_rh
      call read_input_file_initialize_distribution_rh

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_initialize_distribution_rh

         implicit none
         
         kxmax = 1e+100
         kxmin = 0.0
         imfac = 0.0
         refac = 1.0

      end subroutine set_default_parameters_initialize_distribution_rh

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_initialize_distribution_rh

         use file_utils, only: input_unit_exist
         implicit none

         namelist /initialize_distribution_rh/ kxmin, kxmax, imfac, refac
         in_file = input_unit_exist("initialize_distribution_rh", dexist)
         if (dexist) read (unit=in_file, nml=initialize_distribution_rh) 

      end subroutine read_input_file_initialize_distribution_rh

   end subroutine read_namelist_initialize_distribution_rh
   
   !****************************************************************************
   !                              RESTART OPTIONS                              !
   !****************************************************************************
   subroutine read_namelist_restart_options(tstart, scale, restart_file, restart_dir, read_many)

      use mp, only: proc0

      implicit none
       
      real, intent(out) :: tstart, scale
      logical, intent(out) :: read_many
      character(len=300), intent(inout) :: restart_file
      character(len=150), intent(inout) :: restart_dir

      if (.not. proc0) return
      call set_default_parameters_restart_options
      call read_input_file_restart_options

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_restart_options
      
         use file_utils, only: run_name

         implicit none
         
         tstart = 0.0
         scale = 1.0
         restart_file = trim(run_name)//".nc"
         restart_dir = "./"
         read_many = .true.
         
         ! Note that the <read_many> and <save_many> are initialized at the
         ! start of this module so that they are accessible to stella_save.fpp

      end subroutine set_default_parameters_restart_options

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_restart_options

         use file_utils, only: input_unit_exist
         implicit none

         namelist /restart_options/ tstart, scale, restart_file, restart_dir, read_many
         in_file = input_unit_exist("restart_options", dexist)
         if (dexist) read (unit=in_file, nml=restart_options) 

      end subroutine read_input_file_restart_options

   end subroutine read_namelist_restart_options

end module input_file

