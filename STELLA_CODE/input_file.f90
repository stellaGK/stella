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
!   init_distribution (renamed from init_g_knobs)
!   init_distribution_default
!   init_distribution_noise
!   init_distribution_kpar
!   init_distribution_rh
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
   
   private

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

end module input_file

