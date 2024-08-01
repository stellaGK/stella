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

   !> EM
   public :: nfields

   !> Global Logicals
   public :: fields_updated
   !>
   public :: time_field_solve

   private

   !> Logicals
   logical :: fields_updated = .false.
   logical :: fields_initialized = .false.

   !> EM 
   integer :: nfields
   !> 
   real, dimension(2, 5) :: time_field_solve = 0.
   logical :: debug = .false.

   !> TODO-GA: MOVE 
   interface get_dchidy
      module procedure get_dchidy_4d
      module procedure get_dchidy_2d
   end interface get_dchidy

contains

!###############################################################################
!########################### ADVANCE FIELDS VMULO ##############################
!###############################################################################

   subroutine advance_fields(g, phi, apar, bpar, dist, implicit_solve)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use dist_fn_arrays, only: gvmu
      use zgrid, only: nzgrid

      use fields_fluxtube, only: advance_fields_fluxtube
      use fields_full_flux_surface, only: get_fields_ffs
      use fields_arrays, only: time_field_solve

      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: full_flux_surface

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      character(*), intent(in) :: dist

      if (fields_updated) return

      !> Time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

      !> TODO-GA: Add full flux and electromagnetic routines
      nfields = 1

      !> Get phi for Full Flux Surface or Fluxtube
      if (full_flux_surface) then
         if (present(implicit_solve)) then
            if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_ffs_const_in_alpha'
            call get_fields_ffs(g, phi, apar, implicit_solve=.true.)
         else
            if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_ffs'
            call get_fields_ffs(g, phi, apar)
         end if
      else
         if (debug) write (*, *) 'fields::advance_fields_vmulo::get_fields_fluxtube_vmlo'
         call advance_fields_fluxtube(g, phi, apar, bpar, dist)
      end if

      !> Set a flag to indicate that the fields have been updated
      !> this helps avoid unnecessary field solves
      fields_updated = .true.
      !> Time the communications + field solve
      if (proc0) call time_message(.false., time_field_solve(:, 1), ' fields')

   end subroutine advance_fields


!###############################################################################
!################################# UPDATE PHI #################################
!###############################################################################

   subroutine get_phi(phi, dist, skip_fsa)

      use mp, only: proc0, mp_abort, job
      use job_manage, only: time_message
      use parameters_physics, only: radial_variation
      use parameters_numerical, only: ky_solve_radial, ky_solve_real
      use zgrid, only: nzgrid, ntubes
      use geometry, only: dl_over_b
      use arrays_kxky, only: nakx, naky, zonal_mode
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      use species, only: spec, has_electron_species
      use multibox, only: mb_get_phi
      use arrays_fields, only: gamtot, gamtot3
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      logical, optional, intent(in) :: skip_fsa
      real, dimension(:, :, :, :), allocatable :: gamtot_t
      integer :: ia, it, ikx
      complex :: tmp
      logical :: skip_fsa_local
      logical :: has_elec, adia_elec
      logical :: global_quasineutrality, center_cell
      logical :: multibox_mode

      character(*), intent(in) :: dist

      if (debug) write (*, *) 'fields::get_phi'

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ia = 1
      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec &
                  .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      global_quasineutrality = radial_variation .and. ky_solve_radial > 0
      multibox_mode = runtype_option_switch == runtype_multibox
      center_cell = multibox_mode .and. job == 1 .and. .not. ky_solve_real

      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')
      if (dist == 'h') then
         phi = phi / gamtot_h
      else if (dist == 'g' .or. dist == 'gbar') then
         if (global_quasineutrality .and. (center_cell .or. .not. multibox_mode) .and. .not. ky_solve_real) then
            call get_phi_radial(phi)
         else if (global_quasineutrality .and. center_cell .and. ky_solve_real) then
            call mb_get_phi(phi, has_elec, adia_elec)
         else
            ! divide <g> by sum_s (\Gamma_0s-1) Zs^2*e*ns/Ts to get phi
            allocate (gamtot_t(naky, nakx, -nzgrid:nzgrid, ntubes))
            gamtot_t = spread(gamtot, 4, ntubes)
            where (gamtot_t < epsilon(0.0))
               phi = 0.0
            elsewhere
               phi = phi / gamtot_t
            end where
            deallocate (gamtot_t)
         end if
      else
         if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
         call mp_abort('unknown dist option in get_fields. aborting')
         return
      end if

      if (any(gamtot(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')

      ! now handle adiabatic electrons if needed
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'get_phi_adia_elec')
      if (adia_elec .and. zonal_mode(1) .and. .not. skip_fsa_local) then
         if (debug) write (*, *) 'dist_fn::advance_stella::adiabatic_electrons'

         if (dist == 'h') then
            do it = 1, ntubes
               do ikx = 1, nakx
                  tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                  phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * gamtot3_h
               end do
            end do
         else if (dist == 'g' .or. dist == 'gbar') then
            if (global_quasineutrality .and. center_cell .and. ky_solve_real) then
               !this is already taken care of in mb_get_phi
            elseif (global_quasineutrality .and. (center_cell .or. .not. multibox_mode) &
                    .and. .not. ky_solve_real) then
               call add_adiabatic_response_radial(phi)
            else
               do ikx = 1, nakx
                  do it = 1, ntubes
                     tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                     phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * gamtot3(ikx, :)
                  end do
               end do
            end if
         else
            if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
            call mp_abort('unknown dist option in get_fields. aborting')
         end if
      end if
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'get_phi_adia_elec')

   end subroutine get_phi

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
      use fields_arrays, only: phi, apar
      use dist_fn_arrays, only: gnew, gvmu
      use volume_averages, only: volume_average
      use job_manage, only: njobs
      use file_utils, only: runtype_option_switch, runtype_multibox

      implicit none

      real, intent(in) :: target_amplitude
      real :: phi2, rescale

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
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      use parameters_physics, only: include_apar
      use parameters_physics, only: include_bpar
      use parameters_numerical, only: fphi
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: vpa, mu
      use arrays_kxky, only: nakx, aky, naky

      use gyro_averages, only: j0_ffs
      use parameters_physics, only: full_flux_surface

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out) :: dchidy

      integer :: ivmu, iv, is, iky, imu
      complex, dimension(:, :, :, :), allocatable :: field, gyro_tmp

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
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      use parameters_physics, only: include_apar
      use parameters_physics, only: include_bpar
      use parameters_numerical, only: fphi
      use species, only: spec
      use vpamu_grids, only: vpa, mu
      use arrays_kxky, only: nakx, aky, naky

      use gyro_averages, only: j0_ffs
      use parameters_physics, only: full_flux_surface

      implicit none

      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidy

      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: field, gyro_tmp

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
      use gyro_averages, only: gyro_average
      use gyro_averages, only: gyro_average_j1
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: is_idx, iv_idx, imu_idx
      use parameters_physics, only: include_apar
      use parameters_physics, only: include_bpar
      use parameters_numerical, only: fphi
      use species, only: spec
      use vpamu_grids, only: vpa, mu
      use arrays_kxky, only: akx, naky, nakx

      use gyro_averages, only: j0_ffs
      use parameters_physics, only: full_flux_surface

      implicit none

      integer, intent(in) :: ivmu, iz
      complex, dimension(:, :), intent(in) :: phi, apar, bpar
      complex, dimension(:, :), intent(out) :: dchidx

      integer :: iv, is, imu
      complex, dimension(:, :), allocatable :: field, gyro_tmp

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
      use parameters_physics, only: full_flux_surface, radial_variation
      use fields_fluxtube, only: init_fields_fluxtube
      use fields_ffs, only: init_fields_ffs
      use fields_radial_variation, only: init_fields_radial_variation
      implicit none

      debug = debug .and. proc0
      if (fields_initialised) return
      fields_initialised = .true.

      !> Allocate arrays such as phi that are needed throughout the simulation
      if (debug) write (*, *) 'fields::init_fields::allocate_arrays'
      call allocate_arrays

      if (full_flux_surface) then
         if (debug) write (*, *) 'fields::init_fields::init_fields_ffs'
         call init_fields_ffs
      else
         if (debug) write (*, *) 'fields::init_fields::init_fields_fluxtube'
         call init_fields_fluxtube
         nfields = 1

         if (include_apar .or. include_bpar) then 
            if (debug) write (*, *) 'fields::init_fields::init_fields_em'  
            call init_fields_em
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
   !> Allocate arrays needed for solving fields for all versions of stella
   !============================================================================
   
   subroutine allocate_arrays

      use fields_arrays, only: phi, phi_old
      use fields_arrays, only: gamtot, gamtot3
      use zgrid, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo
      use kt_grids, only: naky, nakx
      use species, only: spec, has_electron_species
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg

      use parameters_physics, only: include_apar, include_bpar
      use arrays_fields, only: apar, apar_old
      use arrays_fields, only: bpar, bpar_old

      use arrays_fields, only: gamtot, dgamtotdr, gamtot3
      use arrays_fields, only: gamtot13, gamtot31, gamtot33
      use arrays_fields, only: gamtotinv11, gamtotinv13, gamtotinv31, gamtotinv33

      implicit none

      if (.not. allocated(phi)) then
         allocate (phi(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi = 0.
      end if

      if (.not. allocated(phi_old)) then
         allocate (phi_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi_old = 0.
      end if

      !> TODO-GA: neeed to make this such that it is only for EM stella
      !> Cannot do yet as apar and bpar are integrated into other routines
      !> so need to be allocated. Will try to undo this to save memory and 
      !> CPU time 

      if (.not. allocated(apar)) then
         allocate (apar(naky, nakx, -nzgrid:nzgrid, ntubes))
         apar = 0.
      end if
      if (.not. allocated(apar_old)) then
         allocate (apar_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         apar_old = 0.
      end if
      if (.not. allocated(bpar)) then
         allocate (bpar(naky, nakx, -nzgrid:nzgrid, ntubes))
         bpar = 0.
      end if
      if (.not. allocated(bpar_old)) then
         allocate (bpar_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         bpar_old = 0.
      end if

      if (.not. allocated(gamtot)) then
         allocate (gamtot(naky, nakx, -nzgrid:nzgrid)); gamtot = 0.
      end if

      if (.not. allocated(gamtot3)) then
         if (.not. has_electron_species(spec) &
             .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            allocate (gamtot3(nakx, -nzgrid:nzgrid)); gamtot3 = 0.
         else
            allocate (gamtot3(1, 1)); gamtot3 = 0.
         end if
      end if
            
      if (.not. include_apar) then
         if (.not. allocated(apar_denom)) allocate (apar_denom(1, 1, 1)); apar_denom = 0.
      end if
      if (.not. include_bpar) then
         if (.not. allocated(gamtot33)) allocate (gamtot33(1, 1, 1)); gamtot33 = 0.
         if (.not. allocated(gamtot13)) allocate (gamtot13(1, 1, 1)); gamtot13 = 0.
         if (.not. allocated(gamtot31)) allocate (gamtot31(1, 1, 1)); gamtot31 = 0.
         if (.not. allocated(gamtotinv11)) allocate (gamtotinv11(1, 1, 1)); gamtotinv11 = 0.
         if (.not. allocated(gamtotinv31)) allocate (gamtotinv31(1, 1, 1)); gamtotinv31 = 0.
         if (.not. allocated(gamtotinv13)) allocate (gamtotinv13(1, 1, 1)); gamtotinv13 = 0.
         if (.not. allocated(gamtotinv33))allocate (gamtotinv33(1, 1, 1)); gamtotinv31 = 0.
      end if

   end subroutine allocate_arrays

   !============================================================================
   !============================ FINISH THE FIELDS =============================
   !============================================================================
   subroutine finish_fields

      use fields_arrays, only: phi, phi_old
      use fields_arrays, only: gamtot, gamtot3
      use fields_arrays, only: c_mat, theta
      use parameters_physics, only: full_flux_surface, radial_variation
      use fields_ffs, only: finish_fields_ffs
      use fields_radial_variation, only: finish_radial_fields
      !> TODO-GA: move apar stuff to EM fields
      use fields_arrays, only: apar
      implicit none

      if (allocated(phi)) deallocate (phi)
      if (allocated(phi_old)) deallocate (phi_old)
      if (allocated(gamtot)) deallocate (gamtot)
      if (allocated(gamtot3)) deallocate (gamtot3)

      if (allocated(c_mat)) deallocate (c_mat)
      if (allocated(theta)) deallocate (theta)

      !> TODO-GA: REMOVE
      if (allocated(apar)) deallocate (apar)
      if (full_flux_surface) call finish_fields_ffs
      if (radial_variation) call finish_radial_fields

      !> TODO-GA make EM decoupled from fluxtube -- this will then need an if statement 
      call finish_fields_em
      fields_initialized = .false.

   end subroutine finish_fields

end module fields
