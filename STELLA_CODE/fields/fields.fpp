! Module for advancing and intitialising all fields-related arrays
module fields

   use mpi
   use stella_common_types, only: eigen_type
   use stella_common_types, only: coupled_alpha_type, gam0_ffs_type
   use debug_flags, only: debug => fields_debug

   implicit none

   ! Global Routines
   public :: init_fields, finish_fields

   ! Routines for advancing fields in main routine
   public :: advance_fields

   ! Calculations
   public :: rescale_fields

   ! Global 
   public :: fields_updated
   public :: nfields

   private

   ! Logicals
   logical :: fields_updated = .false.
   logical :: fields_initialized = .false.

   ! EM - to calculate the number of fields that are in the simulation (1 for electrostatic)
   integer :: nfields

   ! For the initialisation 
   logical :: fields_initialised = .false. 

contains

!###############################################################################
!############################## ADVANCE FIELDS #################################
!###############################################################################

   !============================================================================
   !============================== ADVANCE FIELDS ==============================
   !============================================================================
   ! This calls the appropriate routines needed to all fields in the main code.
   ! This routine calls the appropriate update depending on the effects
   ! included in the simulation (e.g. Electrostatic, Full Flux surface effects
   ! or Radiatl Variation effects).
   !============================================================================
   subroutine advance_fields(g, phi, apar, bpar, dist, implicit_solve)

      use mp, only: proc0
      use job_manage, only: time_message
      ! Layouts
      use stella_layouts, only: vmu_lo
      ! Arrays
      use arrays_store_useful, only: time_field_solve
      ! Parameters
      use parameters_physics, only: full_flux_surface
      ! Grids
      use grids_z, only: nzgrid
      ! Routines from other field modules
      use fields_fluxtube, only: advance_fields_fluxtube
      use fields_ffs, only: get_fields_ffs


      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar 
      character(*), intent(in) :: dist
      logical, optional, intent(in) :: implicit_solve
      !-------------------------------------------------------------------------
      if (fields_updated) return

      ! Time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

      ! Do we need Full Flux surface effects?
      if (.not. full_flux_surface) then 
         ! This is the routine for advancing fields in fluxtube.
         ! Note that this will include Electrostatic and Electromagnetic effects
         ! as well as any radial variation effects
         if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_fluxtube'
         call advance_fields_fluxtube(g, phi, apar, bpar, dist)
      else 
         ! This is if Full Flux Surface effects are included
         ! This routine is only needed in the 'implicit_solve' algorithm 
         if (present(implicit_solve)) then
            if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_ffs_const_in_alpha'
            call get_fields_ffs(g, phi, apar, implicit_solve=.true.)
         else
            ! This routine is for advancing the full <phi> field in the code with
            ! FFS effects
            if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_ffs'
            call get_fields_ffs(g, phi, apar)
         end if
      end if

      ! Set a flag to indicate that the fields have been updated
      ! this helps avoid unnecessary field solves
      fields_updated = .true.
      ! Time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

   end subroutine advance_fields


!###############################################################################
!########################## CALCULATIONS FOR FIELDS ############################
!###############################################################################

   !============================================================================
   ! Rescale fields, including the distribution function
   !============================================================================
   subroutine rescale_fields(target_amplitude)

      use mp, only: scope, subprocs, crossdomprocs, sum_allreduce
      use job_manage, only: njobs
      use file_utils, only: runtype_option_switch, runtype_multibox
      ! Arrays
      use arrays_store_fields, only: phi, apar
      use arrays_store_distribution_fn, only: gnew, gvmu
      ! Calculations
      use calculations_volume_averages, only: volume_average

      implicit none

      real, intent(in) :: target_amplitude
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

!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !=========================== INITALIZE THE FIELDS ===========================
   !============================================================================
   subroutine init_fields

      use mp, only: proc0
      use linear_solve, only: lu_decomposition

      ! Parameters
      use parameters_physics, only: full_flux_surface, radial_variation

      ! Routined needed to initialise the different field arrays depending on the 
      ! physics being simulated
      use fields_fluxtube, only: init_fields_fluxtube
      use fields_electromagnetic, only: init_fields_electromagnetic
      use fields_ffs, only: init_fields_ffs
      use fields_radial_variation, only: init_fields_radial_variation

      implicit none
      !-------------------------------------------------------------------------
      debug = debug .and. proc0
      if (fields_initialised) return
      fields_initialised = .true.

      ! Allocate arrays such as phi that are needed throughout the simulation
      if (debug) write (*, *) 'fields::init_fields::allocate_arrays'
      call allocate_arrays

      if (full_flux_surface) then
         if (debug) write (*, *) 'fields::init_fields::init_fields_ffs'
         nfields = 1
         call init_fields_ffs
      else
         if (debug) write (*, *) 'fields::init_fields::init_fields_fluxtube'
         nfields = 1
         call init_fields_fluxtube

         if (debug) write (*, *) 'fields::init_fields::init_fields_electromagnetic'  
         call init_fields_electromagnetic (nfields)

         if (radial_variation) then
            if (debug) write (*, *) 'fields::init_fields::init_fields_radial_variation'
            call init_fields_radial_variation
         end if

      end if

      fields_initialised = .true.

   end subroutine init_fields

   !============================================================================
   !============================= ALLOCATE ARRAYS ==============================
   !============================================================================
   ! Allocate arrays needed for solving fields for all versions of stella
   !============================================================================
   
   subroutine allocate_arrays
       
      ! Parameters
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      use grids_species, only: spec, has_electron_species
      
      ! Arrays to allocate
      use arrays_store_fields, only: phi, phi_old
      use arrays_store_useful, only: gamtot, gamtot3
      
      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      
      ! Time routines
      use arrays_store_useful, only: time_field_solve
      
      implicit none
      
      !-------------------------------------------------------------------------

      ! Time routines
      time_field_solve = 0.
      
      ! Allocate electrostatic arrays on each processor
      if (.not. allocated(phi)) then; allocate (phi(naky, nakx, -nzgrid:nzgrid, ntubes)); phi = 0.; end if
      if (.not. allocated(phi_old)) then; allocate (phi_old(naky, nakx, -nzgrid:nzgrid, ntubes)); phi_old = 0.; end if
      if (.not. allocated(gamtot)) then; allocate (gamtot(naky, nakx, -nzgrid:nzgrid)); gamtot = 0.; end if
      
      ! Allocate <gamtot3> if we have adiabatic field-line-averaged electrons
      if (.not. allocated(gamtot3)) then
         if (.not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            allocate (gamtot3(nakx, -nzgrid:nzgrid)); gamtot3 = 0.
         else
            allocate (gamtot3(1, 1)); gamtot3 = 0.
         end if
      end if

   end subroutine allocate_arrays

   !============================================================================
   !============================ FINISH THE FIELDS =============================
   !============================================================================
   subroutine finish_fields

      ! Parameters
      use parameters_physics, only: full_flux_surface, radial_variation 
      ! Arrays
      use arrays_store_fields, only: phi, phi_old
      use arrays_store_useful, only: gamtot, gamtot3
      ! Routines for deallocating arrays fields depending on the physics being simulated
      use fields_ffs, only: finish_fields_ffs
      use fields_radial_variation, only: finish_radial_fields
      use fields_electromagnetic, only: finish_fields_electromagnetic
      implicit none

      if (allocated(phi)) deallocate (phi)
      if (allocated(phi_old)) deallocate (phi_old)
      if (allocated(gamtot)) deallocate (gamtot)
      if (allocated(gamtot3)) deallocate (gamtot3)

      ! Deallocate EM arrays. Even if not needed these arrays are allocated with size 1
      call finish_fields_electromagnetic
      
      if (full_flux_surface) call finish_fields_ffs
      if (radial_variation) call finish_radial_fields

      fields_initialized = .false.

   end subroutine finish_fields

end module fields
