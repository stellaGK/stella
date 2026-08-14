!###############################################################################
!###############################################################################
!###############################################################################
! 
! This module initialises the terms in the gyrokinetic equation.
! 
!###############################################################################
module gyrokinetic_equation_initialisation

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none
   
   ! Make routines available to other modules
   public :: init_gyrokinetic_equation
   public :: finish_gyrokinetic_equation

   private

   ! Only initialise once
   logical :: initialised_gyrokinetic_equation = .false.

contains

!###############################################################################
!############################ INITIALISE TIME ADVANCE ##########################
!###############################################################################

   subroutine init_gyrokinetic_equation

      ! Physics flags
      use iso_fortran_env, only: output_unit
      use mp, only: proc0
      use dissipation_and_collisions, only: include_collisions
      use parameters_physics, only: radial_variation
      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_physics, only: include_apar      

      ! Initialise the main terms of the gyrokinetic equation
      use calculations_timestep, only: init_cfl
      use gk_drive, only: init_wstar
      use gk_magnetic_drift, only: init_wdrift
      use gk_parallel_streaming, only: init_parallel_streaming
      use gk_mirror, only: init_mirror
      
      ! Initialise optional terms of the gyrokinetic equation
      use dissipation_and_collisions, only: init_collisions
      use gk_flow_shear, only: init_flow_shear
      use gk_radial_variation, only: init_radial_variation
      use gk_sources, only: init_quasineutrality_source
      use gk_sources, only: init_source_timeaverage
      use gk_parallel_nonlinearity, only: init_parallel_nonlinearity
      ! use neoclassical_terms, only: init_neoclassical_terms
      
      ! For NEO's neoclassical corrections. 
      use neoclassical_terms_neo, only: neoclassical_is_enabled
      use gk_neo_chi_terms, only: init_neo_chi_terms    
      use gk_neo_apar_terms, only: init_neo_apar_terms
      use gk_neo_dchidz_terms, only: init_neo_dchidz_terms
      use gk_neo_drive, only: init_wstar1, init_wpol
      use gk_neo_drifts, only: init_neo_curv_drift
                                                                                   
      implicit none

# Sources that don't need preprocessing
set(STELLA_SOURCES_f90
  ../STELLA_CODE/init_stella.f90
  ../STELLA_CODE/arrays/arrays.f90
  ../STELLA_CODE/arrays/arrays_distribution_function.f90
  ../STELLA_CODE/arrays/arrays_fields.f90
  ../STELLA_CODE/arrays/arrays_gyro_averages.f90
  ../STELLA_CODE/arrays/initialise_arrays.f90
  ../STELLA_CODE/arrays/initialise_distribution_function.f90
  ../STELLA_CODE/calculations/calculations_add_explicit_terms.f90
  ../STELLA_CODE/calculations/calculations_checksum.f90
  ../STELLA_CODE/calculations/calculations_finite_differences.f90
  ../STELLA_CODE/calculations/calculations_gyro_averages.f90
  ../STELLA_CODE/calculations/calculations_kxky_derivatives.f90
  ../STELLA_CODE/calculations/calculations_kxky.f90
  ../STELLA_CODE/calculations/calculations_timestep.f90
  ../STELLA_CODE/calculations/calculations_tofrom_ghf.f90
  ../STELLA_CODE/calculations/calculations_transforms.f90
  ../STELLA_CODE/calculations/calculations_velocity_integrals.f90
  ../STELLA_CODE/calculations/calculations_volume_averages.f90
  ../STELLA_CODE/diagnostics/diagnostics_distribution.f90
  ../STELLA_CODE/diagnostics/diagnostics_fluxes_fluxtube.f90
  ../STELLA_CODE/diagnostics/diagnostics_fluxes_fullfluxsurface.f90 
  ../STELLA_CODE/diagnostics/diagnostics_fluxes_radialvariation.f90
  ../STELLA_CODE/diagnostics/diagnostics_fluxes.f90
  ../STELLA_CODE/diagnostics/diagnostics_stresses.f90
  ../STELLA_CODE/diagnostics/diagnostics_moments.f90
  ../STELLA_CODE/diagnostics/diagnostics_omega.f90
  ../STELLA_CODE/diagnostics/diagnostics_potential.f90
  ../STELLA_CODE/diagnostics/diagnostics_write_radial_grid.f90
  ../STELLA_CODE/diagnostics/diagnostics.f90
  ../STELLA_CODE/dissipation/collisions_dougherty.f90
  ../STELLA_CODE/dissipation/collisions_fokkerplanck.f90
  ../STELLA_CODE/dissipation/dissipation_and_collisions.f90
  ../STELLA_CODE/dissipation/dissipation_hyper.f90
  ../STELLA_CODE/geometry/geometry_inputprofiles_interface.f90
  ../STELLA_CODE/geometry/geometry_miller.f90
  ../STELLA_CODE/geometry/geometry_vmec_read_netCDF_file.f90
  ../STELLA_CODE/geometry/geometry_vmec.f90
  ../STELLA_CODE/geometry/geometry_zpinch.f90
  ../STELLA_CODE/geometry/geometry.f90
  ../STELLA_CODE/grids/grids_extended_zgrid.f90
  ../STELLA_CODE/grids/grids_kxky.f90
  ../STELLA_CODE/grids/grids_species.f90
  ../STELLA_CODE/grids/grids_species_from_euterpe.f90
  ../STELLA_CODE/grids/grids_time.f90
  ../STELLA_CODE/grids/grids_velocity.f90
  ../STELLA_CODE/grids/grids_z.f90
  ../STELLA_CODE/field_equations/field_equations.f90
  ../STELLA_CODE/field_equations/field_equations_collisions.f90
  ../STELLA_CODE/field_equations/field_equations_electromagnetic.f90
  ../STELLA_CODE/field_equations/field_equations_fluxtube.f90
  ../STELLA_CODE/field_equations/field_equations_fullfluxsurface.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_drive.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_ffs_solve.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_flow_shear.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_implicit_terms.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_magnetic_drift.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_mirror.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_ExB_nonlinearity.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_parallel_nonlinearity.f90
  ../STELLA_CODE/gyrokinetic_equation/gk_parallel_streaming.f90
  ../STELLA_CODE/gyrokinetic_equation/gyrokinetic_equation_explicit.f90
  ../STELLA_CODE/gyrokinetic_equation/gyrokinetic_equation_implicit.f90
  ../STELLA_CODE/gyrokinetic_equation/gyrokinetic_equation_initialisation.f90
  ../STELLA_CODE/neoclassical/neoclassical_terms.f90
  ../STELLA_CODE/parameters/debug_flags.f90
  ../STELLA_CODE/parameters/parameters_diagnostics.f90
  ../STELLA_CODE/parameters/parameters_multibox.f90
  ../STELLA_CODE/parameters/parameters_numerical.f90
  ../STELLA_CODE/parameters/parameters_physics.f90
  ../STELLA_CODE/parallelisation/redistribute.f90
  ../STELLA_CODE/parallelisation/initialise_redistribute.f90
  ../STELLA_CODE/parallelisation/common_types.f90
  ../STELLA_CODE/parallelisation/parallelisation_layouts.f90
  ../STELLA_CODE/parallelisation/timers.f90
  ../STELLA_CODE/radial_variation/multibox.f90
  ../STELLA_CODE/radial_variation/gk_radial_variation.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_debug.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_diagnostics.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_dissipation.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_flow_shear.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_geometry.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_initialise_distribution_function.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_kxky_grid.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_neoclassical_input.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_numerical_parameters.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_parallelisation.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_physics_parameters.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_radial_variation.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_species.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_velocity_grids.f90
  ../STELLA_CODE/read_namelists_from_input_file/namelist_z_grid.f90
  ../STELLA_CODE/utils/define.inc
)

      ! Only initialise once
      if (initialised_gyrokinetic_equation) return
      initialised_gyrokinetic_equation = .true.
      
      ! Set up neoclassical corrections to the equilibrium Maxwellian. This is only
      ! calculated/needed when simulating higher order terms in rhostar for intrinsic rotation.
      ! This should only be computed if sfincs or NEO is selected as the neoclassical option.
      ! if (debug) write (6, *) 'time_advance::init_time_advance::init_neoclassical_terms'
      !     call init_neoclassical_terms
      ! end if      
 
      ! Calculate the term multiplying dg/dvpa in the mirror term and set up either the
      ! semi-Lagrange machinery or the tridiagonal matrix to be inverted if solving implicitly
      if (debug) write (6, *) 'time_advance::init_time_advance::init_mirror'
      call init_mirror
      
      ! Calculate the term multiplying dg/dz in the parallel streaming term
      ! and set up the tridiagonal matrix to be inverted if solving implicitly
      if (debug) write (6, *) 'time_advance::init_time_advance::init_parstream'
      call init_parallel_streaming
      
      ! Allocate and calculate the factors multiplying dg/dx, dg/dy, dphi/dx 
      ! and dphi/dy in the magnetic drift terms
      if (debug) write (6, *) 'time_advance::init_time_advance::init_wdrift'
      call init_wdrift
      
      ! Allocate and calculate the factor multiplying dphi/dy in the gradient drive term
      if (debug) write (6, *) 'time_advance::init_time_advance::init_wstar'
      call init_wstar

      ! If NEO's neoclassical corrections are enabled, then ...
      if (neoclassical_is_enabled()) then
          ! Allocate and calculate the coeffecient multiplying chi in NEO's neoclassical corrections. 
          call init_neo_chi_terms

          ! If apar is included, allocate and calculate the coeffecient multiplying apar.
          if (include_apar) then          
              call init_neo_apar_terms
          end if 

          ! Allocate and calculate the coeffecient multiplying dchi/dz in NEO's neoclassical corrections.
          call init_neo_dchidz_terms

          ! Allocate and calculate the neoclassical counterpart to wstar and the poloidal gradient drive. 
          call init_wstar1
          call init_wpol

          ! Allocate and calculate the neoclassical curvature drift corrections. 
          call init_neo_curv_drift
      end if
      
      ! Calculate the frequency omega_{zeta,k,s} associated with the parallel flow 
      ! shear and save it as <prl_shear>. Calculate the arrays needed for the discrete 
      ! wavenumber-shift method formulated by Hammett for the perpendicular flow shear.
      ! TODO - check the explanation above.
      if (debug) write (6, *) 'time_advance::init_time_advance::init_flow_shear'
      call init_flow_shear
      
      ! ...
      if (debug) write (6, *) 'time_advance::init_time_advance::init_parallel_nonlinearity'
      if (include_parallel_nonlinearity) call init_parallel_nonlinearity 
      
      ! ...
      if (debug) write (6, *) 'time_advance::init_time_advance::init_radial_variation'
      if (radial_variation) call init_radial_variation
      
      ! ...
      if (include_collisions) then
         if (debug) write (6, *) 'time_advance::init_time_advance::init_collisions'
         call init_collisions
      end if
      
      ! ...
      if (debug) write (6, *) 'time_advance::init_time_advance::init_cfl'
      call init_cfl

      ! ...
      if (debug) write (6, *) 'time_advance::init_time_advance::init_source_timeaverage'
      call init_source_timeaverage
      
      ! ...
      if (debug) write (6, *) 'time_advance::init_time_advance::init_quasineutrality_source'
      call init_quasineutrality_source

   end subroutine init_gyrokinetic_equation

!###############################################################################
!####################### FINISH TIME ADVANCE SUBROUTINE ########################
!###############################################################################

   subroutine finish_gyrokinetic_equation

      use calculations_transforms, only: finish_transforms
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: include_apar
      use grids_extended_zgrid, only: finish_extended_zgrid
      use gk_parallel_streaming, only: finish_parallel_streaming
      use gk_mirror, only: finish_mirror
      use gk_flow_shear, only: finish_flow_shear
      use gk_drive, only: finish_wstar
      use gk_magnetic_drift, only: finish_wdrift
      use gk_parallel_nonlinearity, only: finish_parallel_nonlinearity
      use gk_radial_variation, only: finish_radial_variation 
      use dissipation_and_collisions, only: finish_dissipation
      
      ! For NEO's neoclassical corrections. 

      use neoclassical_terms_neo, only: neoclassical_is_enabled, finish_neoclassical_terms_neo
      use gk_neo_chi_terms, only: finish_neo_chi_terms
      use gk_neo_apar_terms, only: finish_neo_apar_terms
      use gk_neo_dchidz_terms, only: finish_neo_dchidz_terms
      use gk_neo_drive, only: finish_wstar1, finish_wpol
      use gk_neo_drifts, only: finish_neo_curv_drift

      implicit none

      !-------------------------------------------------------------------------

      if (full_flux_surface) call finish_transforms
      call finish_dissipation
      call finish_parallel_nonlinearity
      call finish_wstar
      call finish_wdrift
      call finish_radial_variation
      call finish_parallel_streaming
      call finish_flow_shear
      call finish_mirror
    
      ! If NEO's neoclassical corrections are enabled, then ...
      if (neoclassical_is_enabled()) then
          call finish_neoclassical_terms_neo
          call finish_neo_chi_terms

          ! If apar is included, deallocate the apar neoclassical terms.  
          if (include_apar) then
              call finish_neo_apar_terms
          end if 

          call finish_neo_dchidz_terms
          call finish_wstar1
          call finish_wpol
          call finish_neo_curv_drift
      end if

      initialised_gyrokinetic_equation = .false.

   end subroutine finish_gyrokinetic_equation

end module gyrokinetic_equation_initialisation
