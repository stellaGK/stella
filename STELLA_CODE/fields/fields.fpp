!> Module for advancing and intitialising all fields-related arrays
module fields

   use mpi
   use common_types, only: eigen_type
   use common_types, only: coupled_alpha_type, gam0_ffs_type
   use debug_flags, only: debug => fields_debug

   implicit none

   !> Global Routines
   public :: init_fields, finish_fields

   !> Routines for advancing fields in main routine
   public :: advance_fields

   !> TODO-GA: move to separate routine
   !> Calculations
   public :: rescale_fields
   public :: get_dchidy, get_dchidx

   !> Global 
   public :: fields_updated
   public :: nfields

   private

   !> Logicals
   logical :: fields_updated = .false.
   logical :: fields_initialized = .false.

   !> EM - to calculate the number of fields that are in the simulation (1 for electrostatic)
   integer :: nfields

   !> For the initialisation 
   logical :: fields_initialised = .false. 

   !> TODO-GA: MOVE 
   interface get_dchidy
      module procedure get_dchidy_4d
      module procedure get_dchidy_2d
   end interface get_dchidy

contains

!###############################################################################
!############################## ADVANCE FIELDS #################################
!###############################################################################

   !============================================================================
   !============================== ADVANCE FIELDS ==============================
   !============================================================================
   !> This calls the appropriate routines needed to all fields in the main code.
   !> This routine calls the appropriate update depending on the effects
   !> included in the simulation (e.g. Electrostatic, Full Flux surface effects
   !> or Radiatl Variation effects).
   !============================================================================
   subroutine advance_fields(g, phi, apar, bpar, dist, implicit_solve)

      use mp, only: proc0
      use job_manage, only: time_message
      !> Layouts
      use stella_layouts, only: vmu_lo
      !> Arrays
      use arrays_fields, only: time_field_solve
      !> Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface
      !> Grids
      use zgrid, only: nzgrid
      !> Routines from other field modules
      use fields_fluxtube, only: advance_fields_fluxtube
      use fields_ffs, only: get_fields_ffs


      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar 
      character(*), intent(in) :: dist
      logical, optional, intent(in) :: implicit_solve
      !-------------------------------------------------------------------------
      if (fields_updated) return

      !> Time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

      !> Do we need Full Flux surface effects?
      if (.not. full_flux_surface) then 
         !> This is the routine for advancing fields in fluxtube.
         !> Note that this will include Electrostatic and Electromagnetic effects
         !> as well as any radial variation effects
         if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_fluxtube'
         call advance_fields_fluxtube(g, phi, apar, bpar, dist)
      else 
         !> This is if Full Flux Surface effects are included
         !> This routine is only needed in the 'implicit_solve' algorithm 
         if (present(implicit_solve)) then
            if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_ffs_const_in_alpha'
            call get_fields_ffs(g, phi, apar, implicit_solve=.true.)
         else
            !> This routine is for advancing the full <phi> field in the code with
            !> FFS effects
            if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_ffs'
            call get_fields_ffs(g, phi, apar)
         end if
      end if

      !> Set a flag to indicate that the fields have been updated
      !> this helps avoid unnecessary field solves
      fields_updated = .true.
      !> Time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

   end subroutine advance_fields


!###############################################################################
!########################## CALCULATIONS FOR FIELDS ############################
!###############################################################################

   !============================================================================
   !============================ FIELDS DERIVATIVES ============================
   !============================================================================
   !> Rescale fields, including the distribution function
   !============================================================================
   subroutine rescale_fields(target_amplitude)

      use mp, only: scope, subprocs, crossdomprocs, sum_allreduce
      use job_manage, only: njobs
      use file_utils, only: runtype_option_switch, runtype_multibox
      !> Arrays
      use arrays_fields, only: phi, apar
      use arrays_dist_fn, only: gnew, gvmu
      !> Calculations
      use volume_averages, only: volume_average

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

   !============================================================================
   !============================ FIELDS DERIVATIVES ============================
   !============================================================================
   !> Compute d<chi>/dy and d<chi>/dx in (ky,kx) space where <.> is a gyroaverage
   !>    d<chi>/dy = i * ky * J0 * chi
   !>    d<chi>/dx = i * kx * J0 * chi
   !>    chi = phi - Z/T * vpa * apar
   !> There are different routines depending on the size of the input array
   !============================================================================

   !> TODO-GA: maybe separate for EM and electrostatic
   !> Compute d<chi>/dy in (ky,kx,z,tube) space
   subroutine get_dchidy_4d(phi, apar, bpar, dchidy)

      use constants, only: zi
      !> Layouts
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      !> Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface
      use parameters_numerical, only: fphi
      use parameters_kxky_grids, only: nakx, naky
      !> Grids
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: vpa, mu
      use grids_kxky, only: aky
      !> Calculations
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use gyro_averages, only: j0_ffs

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dchidy

      integer :: ivmu, iv, is, iky, imu
      complex, dimension(:, :, :, :), allocatable :: field, gyro_tmp
      !-------------------------------------------------------------------------
      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (gyro_tmp(naky, nakx, -nzgrid:nzgrid, ntubes))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         is = is_idx(vmu_lo, ivmu)
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         ! intermediate calculation to get factor involving phi contribution
         field = fphi * phi
         ! add apar contribution if including it
         if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
         ! take spectral y-derivative
         do iky = 1, naky
            field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
         end do
         if (full_flux_surface) then
            call gyro_average(field, dchidy(:, :, :, :, ivmu), j0_ffs(:, :, :, ivmu))
         else
            call gyro_average(field, ivmu, dchidy(:, :, :, :, ivmu))
         end if
         if (include_bpar) then
            field = 4.0 * mu(imu) * (spec(is)%tz) * bpar
            do iky = 1, naky
               field(iky, :, :, :) = zi * aky(iky) * field(iky, :, :, :)
            end do
            call gyro_average_j1(field, ivmu, gyro_tmp)
            !> include bpar contribution
            dchidy(:, :, :, :, ivmu) = dchidy(:, :, :, :, ivmu) + gyro_tmp
         end if
      end do

      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidy_4d

   !> Compute d<chi>/dy in (ky,kx) space
   subroutine get_dchidy_2d(iz, ivmu, phi, apar, bpar, dchidy)

      use constants, only: zi
      !> Layouts
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      !> Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface
      use parameters_numerical, only: fphi
      use parameters_kxky_grids, only: nakx, naky
      !> Grids
      use species, only: spec
      use vpamu_grids, only: vpa, mu
      use grids_kxky, only: aky
      !> Calculations
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use gyro_averages, only: j0_ffs

      implicit none

      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidy

      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: field, gyro_tmp
      !-------------------------------------------------------------------------
      allocate (field(naky, nakx))
      allocate (gyro_tmp(naky, nakx))

      is = is_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      field = fphi * phi
      if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
      field = zi * spread(aky, 2, nakx) * field

      if (full_flux_surface) then
         call gyro_average(field, dchidy, j0_ffs(:, :, iz, ivmu))
      else
         call gyro_average(field, iz, ivmu, dchidy)
      end if

      if (include_bpar) then
         field = 4.0 * mu(imu) * (spec(is)%tz) * bpar
         field = zi * spread(aky, 2, nakx) * field
         call gyro_average_j1(field, iz, ivmu, gyro_tmp)
         !> include bpar contribution
         dchidy = dchidy + gyro_tmp
      end if
      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidy_2d

   !> Compute d<chi>/dx in (ky,kx) space
   subroutine get_dchidx(iz, ivmu, phi, apar, bpar, dchidx)

      use constants, only: zi
      !> Layouts
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      !> Parameters
      use parameters_kxky_grids, only: naky, nakx
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface
      use parameters_numerical, only: fphi
      !> Grids
      use species, only: spec
      use vpamu_grids, only: vpa, mu
      use grids_kxky, only: akx
      !> Calculations
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use gyro_averages, only: j0_ffs

      implicit none

      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidx

      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: field, gyro_tmp
      !-------------------------------------------------------------------------
      allocate (field(naky, nakx))
      allocate (gyro_tmp(naky, nakx))

      is = is_idx(vmu_lo, ivmu)
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      field = fphi * phi
      if (include_apar) field = field - 2.0 * vpa(iv) * spec(is)%stm_psi0 * apar
      field = zi * spread(akx, 1, naky) * field

      if (full_flux_surface) then
         call gyro_average(field, dchidx, j0_ffs(:, :, iz, ivmu))
      else
         call gyro_average(field, iz, ivmu, dchidx)
      end if

      if (include_bpar) then
         field = 4 * mu(imu) * (spec(is)%tz) * bpar
         field = zi * spread(akx, 1, naky) * field
         call gyro_average_j1(field, iz, ivmu, gyro_tmp)
         !> include bpar contribution
         dchidx = dchidx + gyro_tmp
      end if
      deallocate (field)
      deallocate (gyro_tmp)

   end subroutine get_dchidx

!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !=========================== INITALIZE THE FIELDS ===========================
   !============================================================================
   subroutine init_fields

      use mp, only: proc0
      use linear_solve, only: lu_decomposition

      !> Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface, radial_variation

      !> Routined needed to initialise the different field arrays depending on the 
      !> physics being simulated
      use fields_fluxtube, only: init_fields_fluxtube
      use fields_electromagnetic, only: init_fields_electromagnetic
      use fields_ffs, only: init_fields_ffs
      use fields_radial_variation, only: init_fields_radial_variation

      implicit none
      !-------------------------------------------------------------------------
      debug = debug .and. proc0
      if (fields_initialised) return
      fields_initialised = .true.

      !> Allocate arrays such as phi that are needed throughout the simulation
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

         if (include_apar .or. include_bpar) then 
            if (debug) write (*, *) 'fields::init_fields::init_fields_electromagnetic'  
            call init_fields_electromagnetic (nfields)
         end if 

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
      
      ! Grids
      use zgrid, only: nzgrid, ntubes
      use parameters_kxky_grids, only: naky, nakx
      
      ! Time routines
      use arrays_fields, only: time_field_solve
      
      implicit none
      
      !-------------------------------------------------------------------------

      ! Time routines
      time_field_solve = 0.
      
      ! Allocate the electrostatic and electromagnetic arrays vs (kx,ky,z,ntubes)
      call allocate_arrays_electrostatic()
      
      ! TODO REMOVE THIS
      call allocate_arrays_electromagnetic()

   contains
      
      !========================== ELECTROSTATIC ARRAYS =========================
      subroutine allocate_arrays_electrostatic
       
         ! Parameters
         use parameters_physics, only: adiabatic_option_switch
         use parameters_physics, only: adiabatic_option_fieldlineavg
         use species, only: spec, has_electron_species
         
         ! Arrays to allocate
         use arrays_fields, only: phi, phi_old
         use arrays_fields, only: gamtot, gamtot3
      
         implicit none
         
         !----------------------------------------------------------------------
         
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
      
      end subroutine allocate_arrays_electrostatic
      
      !========================= ELECTROMAGNETIC ARRAYS ========================
      subroutine allocate_arrays_electromagnetic
      
         ! Parameters 
         use parameters_physics, only: include_apar, include_bpar
      
         ! Arrays to allocate
         use arrays_fields, only: apar, apar_old
         use arrays_fields, only: bpar, bpar_old
         use arrays_fields, only: gamtot13, gamtot31, gamtot33
         use arrays_fields, only: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33
         use arrays_fields, only: apar_denom
      
         implicit none
         
         !----------------------------------------------------------------------

         ! Allocate electromagnetic arrays on each processor
         !if (.not. allocated(apar)) then; allocate (apar(naky, nakx, -nzgrid:nzgrid, ntubes)); apar = 0.; end if
         !if (.not. allocated(bpar)) then; allocate (bpar(naky, nakx, -nzgrid:nzgrid, ntubes)); bpar = 0.; end if
         !if (.not. allocated(apar_old)) then; allocate (apar_old(naky, nakx, -nzgrid:nzgrid, ntubes)); apar_old = 0.; end if
         !if (.not. allocated(bpar_old)) then; allocate (bpar_old(naky, nakx, -nzgrid:nzgrid, ntubes)); bpar_old = 0.; end if
         
         if (.not. include_apar) then
            if (.not. allocated(apar_denom)) then; allocate (apar_denom(1, 1, 1)); apar_denom = 0.; end if
         end if
         
         if (.not. include_bpar) then
            if (.not. allocated(gamtot33)) then; allocate (gamtot33(1, 1, 1)); gamtot33 = 0.; end if
            if (.not. allocated(gamtot13)) then; allocate (gamtot13(1, 1, 1)); gamtot13 = 0.; end if
            if (.not. allocated(gamtot31)) then; allocate (gamtot31(1, 1, 1)); gamtot31 = 0.; end if
            if (.not. allocated(gamtotinv11)) then; allocate (gamtotinv11(1, 1, 1)); gamtotinv11 = 0.; end if
            if (.not. allocated(gamtotinv31)) then; allocate (gamtotinv31(1, 1, 1)); gamtotinv31 = 0.; end if
            if (.not. allocated(gamtotinv13)) then; allocate (gamtotinv13(1, 1, 1)); gamtotinv13 = 0.; end if
            if (.not. allocated(gamtotinv33)) then; allocate (gamtotinv33(1, 1, 1)); gamtotinv31 = 0.; end if
         end if
      
      end subroutine allocate_arrays_electromagnetic

   end subroutine allocate_arrays

   !============================================================================
   !============================ FINISH THE FIELDS =============================
   !============================================================================
   subroutine finish_fields

      !> Parameters
      use parameters_physics, only: full_flux_surface, radial_variation 
      use parameters_physics, only: include_apar, include_bpar    
      !> Arrays
      use arrays_fields, only: phi, phi_old
      use arrays_fields, only: gamtot, gamtot3
      !> TODO-GA: move apar stuff to EM fields
      use arrays_fields, only: apar, apar_denom
      use arrays_fields, only: apar_old, bpar_old
      use arrays_fields, only: gamtot, gamtot3
      use arrays_fields, only: gamtot13, gamtot31, gamtot33
      use arrays_fields, only: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33
      use arrays_fields, only: apar_denom
      !> Routines for deallocating arrays fields depending on the physics being simulated
      use fields_ffs, only: finish_fields_ffs
      use fields_radial_variation, only: finish_radial_fields
      use fields_electromagnetic, only: finish_fields_electromagnetic
      implicit none

      if (allocated(phi)) deallocate (phi)
      if (allocated(phi_old)) deallocate (phi_old)
      if (allocated(gamtot)) deallocate (gamtot)
      if (allocated(gamtot3)) deallocate (gamtot3)

      !> TODO-GA: REMOVE
      if (allocated(apar_old)) deallocate(apar_old)
      if (allocated(bpar_old)) deallocate(bpar_old)
      if (allocated(apar)) deallocate (apar)
      if (allocated(apar_denom)) deallocate (apar_denom)
      if (allocated(gamtot33)) deallocate (gamtot33)
      if (allocated(gamtot13)) deallocate (gamtot13)
      if (allocated(gamtot31)) deallocate (gamtot31)
      if (allocated(gamtotinv11)) deallocate(gamtotinv11)
      if (allocated(gamtotinv31)) deallocate(gamtotinv31)
      if (allocated(gamtotinv13)) deallocate(gamtotinv13)
      if (allocated(gamtotinv33)) deallocate(gamtotinv33)

      !> TODO-GA: move the above deallocations into 'finish_fields_electromagnetic' when 
      !> EM is decoupled 
      if (include_apar .or. include_bpar) call finish_fields_electromagnetic
      if (full_flux_surface) call finish_fields_ffs
      if (radial_variation) call finish_radial_fields

      fields_initialized = .false.

   end subroutine finish_fields

end module fields
