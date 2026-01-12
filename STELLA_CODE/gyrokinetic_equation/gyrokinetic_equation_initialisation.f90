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
      use gk_drive, only: init_wstar, init_wpol
      use gk_magnetic_drift, only: init_wdrift
      use gk_parallel_streaming, only: init_parallel_streaming
      use gk_mirror, only: init_mirror
      
      ! Initialise optional terms of the gyrokinetic equation
      use dissipation_and_collisions, only: init_collisions
      use gk_flow_shear, only: init_flow_shear
      use gk_radial_variation, only: init_radial_variation
      use gk_sources, only: init_quasineutrality_source
      use gk_sources, only: init_source_timeaverage
      use gk_nonlinearity, only: init_parallel_nonlinearity
      ! use neoclassical_terms, only: init_neoclassical_terms
      use neoclassical_terms_neo, only: neoclassical_is_enabled, init_neoclassical_terms_neo
      use gk_neoclassical_chi_terms, only: init_neoclassical_chi_terms    
      use gk_neoclassical_apar_terms, only: init_neoclassical_apar_terms
                                                                                   
      implicit none

      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_gyrokinetic_equation) return
      initialised_gyrokinetic_equation = .true.
      
      ! Set up neoclassical corrections to the equilibrium Maxwellian. This is only
      ! calculated/needed when simulating higher order terms in rhostar for intrinsic rotation.
      ! This should only be computed if sfincs or NEO is selected as the neoclassical option.
      ! if (debug) write (6, *) 'time_advance::init_time_advance::init_neoclassical_terms'
      !     call init_neoclassical_terms
      ! end if      

      ! Set up neoclassical corrections to the equilibrium Maxwellian. This is only
      ! calculated/needed when simulating higher order terms in rhostar for intrinsic rotation or steep gradient regimes.                          
      ! This option should only be computed if NEO is selected as the neoclassical option (neo_option_switch = 2).
      if (debug) write (6, *) 'time_advance::init_time_advance::init_neoclassical_terms'
      if (neoclassical_is_enabled()) then
          call init_neoclassical_terms_neo                                        
          if (proc0) then
              write(output_unit, '(A)') '!==================================!'
              write(output_unit, '(A)') '  NEOs Neoclassical terms enabled.  '          
              write(output_unit, '(A)') '!==================================!' 
          end if
      end if 
 
     
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
          ! NOT CURRENTLY WORKING.
          ! call init_wpol

          ! Allocate and calculate the coeffecient multiplying phi in NEO's neoclassical corrections. 
          call init_neoclassical_chi_terms

          ! Allocate and calculate the coeffecient multiplying dphi/dz in NEO's neoclassical corrections.
          ! call init_neoclassical_dchi/dz_terms

          ! If apar is included, allocate and calculate the coeffecient multiplying apar.
          if (include_apar) then          
              call init_neoclassical_apar_terms
          end if 
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
      use gk_drive, only: finish_wstar, finish_wpol
      use gk_magnetic_drift, only: finish_wdrift
      use gk_nonlinearity, only: finish_parallel_nonlinearity
      use gk_radial_variation, only: finish_radial_variation 
      use dissipation_and_collisions, only: finish_dissipation
      
      ! For NEO's neoclassical corrections. 

      use neoclassical_terms_neo, only: neoclassical_is_enabled, finish_neoclassical_terms_neo
      use gk_neoclassical_chi_terms, only: finish_neoclassical_chi_terms
      use gk_neoclassical_apar_terms, only: finish_neoclassical_apar_terms

      implicit none

      !-------------------------------------------------------------------------

      if (full_flux_surface) call finish_transforms
      call finish_dissipation
      call finish_parallel_nonlinearity
      call finish_wstar
      call finish_wpol
      call finish_wdrift
      call finish_radial_variation
      call finish_parallel_streaming
      call finish_flow_shear
      call finish_mirror
    
      ! If NEO's neoclassical corrections are enabled, then ...

      if (neoclassical_is_enabled()) then
          call finish_neoclassical_terms_neo

          ! NOT CURRENTLY WORKING.
          ! finish init_wpol
          call finish_neoclassical_chi_terms
          ! call finish_neoclassical_dchi/dz_terms

          ! If apar is included, deallocate the apar neoclassical terms.  

          if (include_apar) then
              call finish_neoclassical_apar_terms
          end if 

          ! call finish_neoclassical_apar_terms
      end if

      initialised_gyrokinetic_equation = .false.

   end subroutine finish_gyrokinetic_equation

end module gyrokinetic_equation_initialisation
