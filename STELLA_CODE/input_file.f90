!###############################################################################
!########################## READ ALL STELLA NAMELISTS ##########################
!###############################################################################
! 
! This module will set the default input input parameters for each name list,
! and it will read the stella input file per namelist.
! 
! For each namelists two routines will exist:
!    - read_namelist_<namelist>
!    - set_default_parameters_<namelist>
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
!   initialize_potential (renamed from init_g_knobs)
!   initialize_potential_default
!   initialize_potential_noise
!   initialize_potential_kpar
!   initialize_potential_rh
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
   public :: read_namelist_initialize_potential
   
   ! Parameters need to be public
   public :: init_potential_option_default, init_potential_option_noise, init_potential_option_restart_many
   public :: init_potential_option_kpar, init_potential_option_rh, init_potential_option_remap
   
   private
   
   ! Create parameters for the <init_potential_option>
   integer, parameter :: init_potential_option_default = 1
   integer, parameter :: init_potential_option_noise = 2
   integer, parameter :: init_potential_option_restart_many = 3
   integer, parameter :: init_potential_option_kpar = 4
   integer, parameter :: init_potential_option_rh = 5
   integer, parameter :: init_potential_option_remap = 6

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
      call set_default_parameters_dissipation
      call read_input_file_dissipation
      call check_inputs_dissipation

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_dissipation

         implicit none

         ! By default we do not include collisions nor dissipation
         include_collisions = .false.
         collisions_implicit = .true.
         hyper_dissipation = .false.

         ! Options: dougherty or fokker-planck
         collision_model = "dougherty"

      end subroutine set_default_parameters_dissipation

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_dissipation

         use file_utils, only: input_unit_exist

         implicit none

         integer :: in_file
         logical :: dexist

         ! Variables in the <dissipation> namelist
         namelist /dissipation/ include_collisions, collisions_implicit, collision_model, hyper_dissipation

         ! Overwrite the default input parameters by those specified in the input file
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
   !                           INITIALIZE POTENTIAL                            !
   !****************************************************************************
   subroutine read_namelist_initialize_potential(init_potential_switch, phiinit, left, chop_side, scale_to_phiinit)

      use mp, only: proc0

      implicit none
      
      ! Variables that are read from the input file
      real, intent(out) :: phiinit
      logical, intent(out) :: left
      logical, intent(out) :: chop_side
      logical, intent(out) :: scale_to_phiinit
      integer, intent(out) :: init_potential_switch
      
      ! Local variable to set <init_potential_switch>
      character(20) :: initialize_potential_option

      if (.not. proc0) return
      call set_default_parameters_initialize_potential
      call read_input_file_initialize_potential

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_initialize_potential

         implicit none

         ! Options: {default, noise, many, kpar, rh, remap}
         initialize_potential_option = "default"
         
         ! Other options
         phiinit = 1.0
         scale_to_phiinit = .false.
         chop_side = .false.
         left = .true.

      end subroutine set_default_parameters_initialize_potential

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_initialize_potential

         use file_utils, only: input_unit_exist, error_unit
         use text_options, only: text_option, get_option_value

         implicit none

         ! Variables needed to read the input file
         integer :: ierr
         integer :: in_file
         logical :: dexist
      
         ! Link text options for <initialize_potential_option> to an integer value
         type(text_option), dimension(6), parameter :: init_potential_options = &
             (/text_option('default', init_potential_option_default), &
               text_option('noise', init_potential_option_noise), &
               text_option('many', init_potential_option_restart_many), &
               text_option('kpar', init_potential_option_kpar), &
               text_option('rh', init_potential_option_rh), &
               text_option('remap', init_potential_option_remap)/)

         ! Variables in the <initialize_potential> namelist
         namelist /initialize_potential/ initialize_potential_option, phiinit, scale_to_phiinit, chop_side, left
         
         !----------------------------------------------------------------------

         ! Overwrite the default input parameters by those specified in the input file
         in_file = input_unit_exist("initialize_potential", dexist)
         if (dexist) read (unit=in_file, nml=initialize_potential)
         
         ! Read the text option in <initialize_potential> and store it in <adiabatic_option_switch>
         ierr = error_unit()
         call get_option_value(initialize_potential_option, init_potential_options, init_potential_switch, &
            ierr, "initialize_potential_option in initialize_potential")

      end subroutine read_input_file_initialize_potential

   end subroutine read_namelist_initialize_potential

end module input_file

