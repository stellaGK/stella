!###############################################################################
!################################ ADVANCE FIELDS ###############################
!###############################################################################
! 
! This module is used to advance the fields. It calls the appropriate routines
! depending on the physics included, and the simulation domain.
! 
! The following fields are evolved within stella:
!   - The perturbed electrostatic potential (phi)
!   - The perturbed parallel magnetic potential (apar)
!   - The perturbed parallel magnetic field (bpar)
! 
! Physics options:
!        - Electrostatic
!        - Electromagnetic
! 
! Simulation domains:
!        - Flux tube -> Supports electrostatic and electromagnetic
!        - Full flux annulus -> Supports only electrostatic
!        - Radially global (extension of flux tube) -> Supports only electrostatic
! 
! The corresponding equations that are solved can be found in the appropriate
! module. 
!###############################################################################
module field_equations

   ! Load the debug flags
   use debug_flags, only: debug => fields_debug

   implicit none

   ! Make routines available to other modules
   public :: init_field_equations
   public :: finish_field_equations
   public :: advance_fields
   public :: rescale_fields

   ! Make parameters available to other modules
   public :: fields_updated
   public :: nfields

   private

   ! Number of fields that are in the simulation
   ! 1 for electrostatic, i.e. phi
   ! 2 or 3 for electromagnetic, i.e., apar and bpar
   integer :: nfields

   ! Every time the distribution function is updated
   ! we want to update the fields as well
   logical :: fields_updated = .false.
   
   ! Only initialise once
   logical :: initialised_fields = .false.

contains

!###############################################################################
!############################## ADVANCE ALL FIELDS #############################
!###############################################################################

   !============================================================================
   !============== ADVANCE ALL FIELDS USING THE FIELD EQUATIONS ================
   !============================================================================
   ! This calls the appropriate routines needed to all fields in the main code.
   ! This routine calls the appropriate update depending on the effects
   ! included in the simulation (e.g. Electrostatic, Full Flux surface effects
   ! or Radiatl Variation effects).
   !============================================================================
   subroutine advance_fields(g, phi, apar, bpar, dist, implicit_solve)

      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo
      use parameters_physics, only: full_flux_surface
      use grids_z, only: nzgrid
      use timers, only: time_field_solve
      
      ! Routines from other field modules
      use field_equations_fluxtube, only: advance_fields_fluxtube
      use field_equations_fullfluxsurface, only: advance_fields_using_field_equations_fullfluxsurface

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      character(*), intent(in) :: dist
      logical, optional, intent(in) :: implicit_solve
      
      !-------------------------------------------------------------------------
      
      ! Only update the fields when it is needed
      ! This flag will be set to False anytime the distribution function is updated
      if (fields_updated) return

      ! Time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

      !-------------------------------------------------------------------------
      !                    Fluxtube simulation (+ Radial Variation)
      !-------------------------------------------------------------------------
      ! Advance the fields for a flux-tube simulation
      ! This will include electrostatic and electromagnetic effects, as well as
      ! any radial variation effects if included.
      !-------------------------------------------------------------------------
      if (.not. full_flux_surface) then 
         if (debug) write (*, *) 'field_equations::advance_fields_fluxtube'
         call advance_fields_fluxtube(g, phi, apar, bpar, dist)

      !-------------------------------------------------------------------------
      !                        Full Flux Surface simulation
      !-------------------------------------------------------------------------
      ! The below routines are only for full-flux-surface simulations
      !-------------------------------------------------------------------------
      else 
      
         ! This routine is only needed in the 'implicit_solve' algorithm
         if (present(implicit_solve)) then
            if (debug) write (*, *) 'field_equations::advance_fields_using_field_equations_fullfluxsurface::implicit'
            call advance_fields_using_field_equations_fullfluxsurface(g, phi, apar, implicit_solve=.true.)
            
         ! This routine is for advancing the full <phi> field in the code with FFS effects
         else
            if (debug) write (*, *) 'field_equations::advance_fields_using_field_equations_fullfluxsurface'
            call advance_fields_using_field_equations_fullfluxsurface(g, phi, apar)
         end if
         
      end if
      !-------------------------------------------------------------------------
      
      ! Set a flag to indicate that the fields have been updated
      ! This helps avoid unnecessary field solves
      fields_updated = .true.
      
      ! Time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

   end subroutine advance_fields

!###############################################################################
!############################ INITALISE & FINALISE #############################
!###############################################################################

   !============================================================================
   !=========================== INITALISE THE FIELDS ===========================
   !============================================================================
   subroutine init_field_equations

      use linear_solve, only: lu_decomposition
      use parameters_physics, only: full_flux_surface, radial_variation

      ! Routines needed to initialise the different field arrays depending on the physics being simulated
      use field_equations_radialvariation, only: init_field_equations_radialvariation
      use field_equations_electromagnetic, only: init_field_equations_electromagnetic
      use field_equations_fluxtube, only: init_field_equations_fluxtube
      use field_equations_fullfluxsurface, only: init_field_equations_fullfluxsurface

      implicit none
      
      !-------------------------------------------------------------------------

      ! Only initialise once
      if (initialised_fields) return
      initialised_fields = .true.

      ! Allocate arrays such as phi that are needed throughout the simulation
      if (debug) write (*, *) 'field_equations::init::allocate_arrays'
      call allocate_arrays

      !-------------------------------------------------------------------------
      !                        Full Flux Surface simulation
      !-------------------------------------------------------------------------
      if (full_flux_surface) then
      
         if (debug) write (*, *) 'field_equations::init::ffs'
         nfields = 1
         call init_field_equations_fullfluxsurface
         
      !-------------------------------------------------------------------------
      !                            Flux Tube simulation
      !-------------------------------------------------------------------------
      else
      
         ! Electrostatic flux tube 
         if (debug) write (*, *) 'field_equations::init::fluxtube'
         nfields = 1
         call init_field_equations_fluxtube

         ! Electromagnetic flux tube
         if (debug) write (*, *) 'field_equations::init::electromagnetic'
         call init_field_equations_electromagnetic(nfields)

         ! Radial Variation effects
         if (radial_variation) then
            if (debug) write (*, *) 'field_equations::init::radial_variation'
            call init_field_equations_radialvariation
         end if
         
      end if

   end subroutine init_field_equations

   !============================================================================
   !============================= ALLOCATE ARRAYS ==============================
   !============================================================================
   ! Allocate arrays needed for solving fields for all versions of stella
   !============================================================================
   subroutine allocate_arrays

      use arrays_fields, only: phi, phi_old
      use arrays, only: denominator_fields, denominator_fields_MBR
      use grids_species, only: adiabatic_option_switch
      use grids_species, only: adiabatic_option_fieldlineavg
      use grids_species, only: spec, has_electron_species
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      
      implicit none
      
      !-------------------------------------------------------------------------
      
      ! Allocate electrostatic arrays on each processor
      if (.not. allocated(phi)) then; allocate (phi(naky, nakx, -nzgrid:nzgrid, ntubes)); phi = 0.; end if
      if (.not. allocated(phi_old)) then; allocate (phi_old(naky, nakx, -nzgrid:nzgrid, ntubes)); phi_old = 0.; end if
      if (.not. allocated(denominator_fields)) then; allocate (denominator_fields(naky, nakx, -nzgrid:nzgrid)); denominator_fields = 0.; end if
      
      ! Allocate <denominator_fields_MBR> if we have adiabatic field-line-averaged electrons
      if (.not. allocated(denominator_fields_MBR)) then
         if (.not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            allocate (denominator_fields_MBR(nakx, -nzgrid:nzgrid)); denominator_fields_MBR = 0.
         else
            allocate (denominator_fields_MBR(1, 1)); denominator_fields_MBR = 0.
         end if
      end if

   end subroutine allocate_arrays

   !============================================================================
   !============================ FINISH THE FIELDS =============================
   !============================================================================
   subroutine finish_field_equations

      use parameters_physics, only: full_flux_surface, radial_variation
      use arrays_fields, only: phi, phi_old
      use arrays, only: denominator_fields, denominator_fields_MBR
      
      ! Routines for deallocating arrays fields depending on the physics being simulated
      use field_equations_fullfluxsurface, only: finish_field_equations_fullfluxsurface
      use field_equations_radialvariation, only: finish_field_equations_radialvariation
      use field_equations_electromagnetic, only: finish_field_equations_electromagnetic
      
      implicit none
      
      !-------------------------------------------------------------------------

      ! Deallocate ararys
      if (allocated(phi)) deallocate (phi)
      if (allocated(phi_old)) deallocate (phi_old)
      if (allocated(denominator_fields)) deallocate (denominator_fields)
      if (allocated(denominator_fields_MBR)) deallocate (denominator_fields_MBR)

      ! Deallocate arrays from other field routines
      call finish_field_equations_electromagnetic
      if (full_flux_surface) call finish_field_equations_fullfluxsurface
      if (radial_variation) call finish_field_equations_radialvariation

      ! The fields are no longer initialised
      initialised_fields = .false.

   end subroutine finish_field_equations
   

!###############################################################################
!########################## CALCULATIONS FOR FIELDS ############################
!###############################################################################

   !============================================================================
   ! Rescale fields, including the distribution function
   !============================================================================
   subroutine rescale_fields(target_amplitude)
  
      use mp, only: scope, subprocs
      use mp, only: crossdomprocs, sum_allreduce
      use job_manage, only: njobs
      use file_utils, only: runtype_option_switch, runtype_multibox
      use arrays_fields, only: phi, apar
      use arrays_distribution_function, only: gnew, gvmu
      use calculations_volume_averages, only: volume_average

      implicit none

      ! Arguments
      real, intent(in) :: target_amplitude
      
      ! Local variables
      real :: phi2, rescale
      
      !-------------------------------------------------------------------------
      
      call volume_average(phi, phi2)

      if (runtype_option_switch == runtype_multibox) then
         call scope(crossdomprocs)
         call sum_allreduce(phi2)
         call scope(subprocs)
         phi2 = phi2 / njobs
      end if

      rescale = target_amplitude / sqrt(phi2)

      phi = rescale * phi
      apar = rescale * apar
      gnew = rescale * gnew
      gvmu = rescale * gvmu

   end subroutine rescale_fields

end module field_equations
