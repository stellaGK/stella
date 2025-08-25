!> Module for advancing and initialising the fields when Radial Variation effects are included
module fields_radial_variation

   use mpi

   implicit none

   !> Advance EM fields routines
   public :: get_radial_correction
   public :: get_phi_for_radial, add_adiabatic_response_radial
   public :: add_radial_correction_int_species
   !> Initialise and Finalise Routines
   !> TODO-GA: probably can make private -- need to do!
   public :: init_fields_radial_variation
   public :: finish_radial_fields

   private

#ifdef ISO_C_BINDING
   integer :: phi_shared_window = MPI_WIN_NULL
#endif

   integer :: zm
   logical :: debug = .false.

contains

!###############################################################################
!####################### ADVANCE RADIALLY GLOBAL FIELDS ########################
!###############################################################################
   subroutine get_phi_for_radial(phi, dist, skip_fsa)

      use mp, only: proc0, mp_abort, job
      use job_manage, only: time_message
      use parameters_physics, only: radial_variation
      use parameters_multibox, only: ky_solve_radial, ky_solve_real
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: zonal_mode
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      use grids_species, only: spec, has_electron_species
      use multibox, only: mb_get_phi
      use arrays_store_fields, only: gamtot
      use file_utils, only: runtype_option_switch, runtype_multibox
      use arrays_store_fields, only: time_field_solve

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      logical, optional, intent(in) :: skip_fsa
      integer :: ia
      logical :: skip_fsa_local
      logical :: has_elec, adia_elec
      logical :: global_quasineutrality, center_cell
      logical :: multibox_mode

      character(*), intent(in) :: dist

      if (debug) write (*, *) 'dist_fn::advance_stella::get_phi'

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
      if (dist == 'gbar') then
         if (global_quasineutrality .and. (center_cell .or. .not. multibox_mode) .and. .not. ky_solve_real) then
            call get_phi_radial(phi)
         else if (global_quasineutrality .and. center_cell .and. ky_solve_real) then
            call mb_get_phi(phi, has_elec, adia_elec)
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

         if (dist == 'gbar') then
            if (global_quasineutrality .and. center_cell .and. ky_solve_real) then
               !this is already taken care of in mb_get_phi
            elseif (global_quasineutrality .and. (center_cell .or. .not. multibox_mode) &
                    .and. .not. ky_solve_real) then
               call add_adiabatic_response_radial(phi)
            end if
         else
            if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
            call mp_abort('unknown dist option in get_fields. aborting')
         end if
      end if
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'get_phi_adia_elec')

   end subroutine get_phi_for_radial

   !> Non-perturbative approach to solving quasineutrality for radially
   !> global simulations
   subroutine get_phi_radial(phi)

#ifdef ISO_C_BINDING
      use mpi
      use mp, only: curr_focus, sharedsubprocs, scope
      use mp, only: split_n_tasks, sgproc0
      use grids_z, only: nztot
      use arrays_store_fields, only: phi_shared
      use mp_lu_decomposition, only: lu_matrix_multiply_local
#endif
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      use parameters_multibox, only: ky_solve_radial
      use grids_z, only: nzgrid, ntubes
      use grids_species, only: spec, has_electron_species
      use parameters_kxky_grid, only: nakx, naky
      use grids_kxky, only: zonal_mode
      use linear_solve, only: lu_back_substitution
      use arrays_store_fields, only: gamtot, phi_solve

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      integer :: it, iz, iky, zmi
      integer :: naky_r
      complex, dimension(:, :), allocatable :: g0k, g0x
      logical :: has_elec, adia_elec
#ifdef ISO_C_BINDING
      integer :: counter, c_lo, c_hi
      integer :: prior_focus, ierr
#endif

      allocate (g0k(1, nakx))
      allocate (g0x(1, nakx))

      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec &
                  .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      naky_r = min(naky, ky_solve_radial)
#ifdef ISO_C_BINDING
      prior_focus = curr_focus
      call scope(sharedsubprocs)

      call split_n_tasks(nztot * ntubes * naky_r, c_lo, c_hi)

      call scope(prior_focus)
      counter = 0
      if (sgproc0) phi_shared = phi
      call mpi_win_fence(0, phi_shared_window, ierr)
#endif
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do iky = 1, naky_r
#ifdef ISO_C_BINDING
               counter = counter + 1
               if ((counter >= c_lo) .and. (counter <= c_hi)) then
                  if (.not. (adia_elec .and. zonal_mode(iky))) then
                     zmi = 0
                     if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
                     call lu_back_substitution(phi_solve(iky, iz)%zloc, &
                                               phi_solve(iky, iz)%idx, phi_shared(iky, (1 + zmi):, iz, it))
                     if (zmi > 0) phi(iky, zmi, iz, it) = 0.0
                  end if
               end if
#else
               if (.not. (adia_elec .and. zonal_mode(iky))) then
                  zmi = 0
                  if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
                  call lu_back_substitution(phi_solve(iky, iz)%zloc, &
                                            phi_solve(iky, iz)%idx, phi(iky, (1 + zmi):, iz, it))
                  if (zmi > 0) phi(iky, zmi, iz, it) = 0.0
               end if
#endif
            end do
         end do
      end do
#ifdef ISO_C_BINDING
      call mpi_win_fence(0, phi_shared_window, ierr)
      phi = phi_shared
#endif

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do iky = naky_r + 1, naky
               phi(iky, :, iz, it) = phi(iky, :, iz, it) / gamtot(iky, :, iz)
            end do
         end do
      end do

      if (ky_solve_radial == 0 .and. any(gamtot(1, 1, :) < epsilon(0.))) &
         phi(1, 1, :, :) = 0.0

      deallocate (g0k, g0x)

   end subroutine get_phi_radial

   !> Add the adiabatic eletron contribution for globally radial simulations.
   !> This actually entails solving for the whole ky = 0 slice of phi at once (not really adding!)
   subroutine add_adiabatic_response_radial(phi)

#ifdef ISO_C_BINDING
      use mpi
      use mp, only: sgproc0, comm_sgroup
      use arrays_store_fields, only: qn_zf_window
      use mp_lu_decomposition, only: lu_matrix_multiply_local
#else
      use linear_solve, only: lu_back_substitution
#endif
      use grids_z, only: nzgrid, ntubes
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use geometry, only: dl_over_b, d_dl_over_b_drho
      use parameters_kxky_grid, only: nakx 
      use parameters_multibox, only: boundary_size
      use grids_kxky, only: rho_d_clamped
      use arrays_store_fields, only: phizf_solve, phi_ext
      use arrays_store_fields, only: phi_proj, phi_proj_stage, theta
      use arrays_store_fields, only: exclude_boundary_regions_qn, exp_fac_qn, tcorr_source_qn

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      integer :: ia, it, iz, ikx
      integer :: inmat
      complex, dimension(:, :), allocatable :: g0k, g1k, g0x
#ifdef ISO_C_BINDING
      integer :: ierr
#endif

      allocate (g0k(1, nakx))
      allocate (g1k(1, nakx))
      allocate (g0x(1, nakx))

      ia = 1

      do it = 1, ntubes
         ! calculate <<g>_psi>_T
         g1k = 0.0
         do iz = -nzgrid, nzgrid - 1
            g0k(1, :) = phi(1, :, iz, it)
            call transform_kx2x_unpadded(g0k, g0x)
            g0x(1, :) = (dl_over_b(ia, iz) + d_dl_over_b_drho(ia, iz) * rho_d_clamped) * g0x(1, :)
            if (exclude_boundary_regions_qn) then
               g0x(1, :) = sum(g0x(1, (boundary_size + 1):(nakx - boundary_size))) &
                           / (nakx - 2 * boundary_size)
               g0x(1, 1:boundary_size) = 0.0
               g0x(1, (nakx - boundary_size + 1):) = 0.0
            else
               g0x(1, :) = sum(g0x(1, :)) / nakx
            end if

            call transform_x2kx_unpadded(g0x, g0k)

            g1k = g1k + g0k
         end do

         phi_proj_stage(:, 1, it) = g1k(1, :)
         if (tcorr_source_qn < epsilon(0.0)) then
            do iz = -nzgrid, nzgrid - 1
               phi(1, :, iz, it) = phi(1, :, iz, it) - g1k(1, :)
            end do
         else
            do iz = -nzgrid, nzgrid - 1
               phi(1, :, iz, it) = phi(1, :, iz, it) &
                                   - (1.-exp_fac_qn) * g1k(1, :) - exp_fac_qn * phi_proj(:, 1, it)
            end do
         end if

#ifdef ISO_C_BINDING
         if (sgproc0) then
#endif
            do iz = -nzgrid, nzgrid - 1
               do ikx = 1, nakx
                  inmat = ikx + nakx * (iz + nzgrid)
                  phi_ext(inmat) = phi(1, ikx, iz, it)
               end do
            end do
#ifdef ISO_C_BINDING
         end if
         call mpi_win_fence(0, qn_zf_window, ierr)
#endif

#ifdef ISO_C_BINDING
         call lu_matrix_multiply_local(comm_sgroup, qn_zf_window, phizf_solve%zloc, phi_ext)
         call mpi_win_fence(0, qn_zf_window, ierr)
#else
         call lu_back_substitution(phizf_solve%zloc, phizf_solve%idx, phi_ext)
#endif

         do iz = -nzgrid, nzgrid - 1
            do ikx = 1, nakx
               inmat = ikx + nakx * (iz + nzgrid)
               phi(1, ikx, iz, it) = phi_ext(inmat)
            end do
         end do

         !enforce periodicity
         phi(1, :, nzgrid, it) = phi(1, :, -nzgrid, it)

         ! calculate Theta.phi
         g1k = 0.0
         do iz = -nzgrid, nzgrid - 1
            do ikx = 1, nakx
               g0k(1, ikx) = sum(theta(ikx, :, iz) * phi(1, :, iz, it))
            end do

            call transform_kx2x_unpadded(g0k, g0x)

            g0x(1, :) = (dl_over_b(ia, iz) + d_dl_over_b_drho(ia, iz) * rho_d_clamped) * g0x(1, :)
            if (exclude_boundary_regions_qn) then
               g0x(1, :) = sum(g0x(1, (boundary_size + 1):(nakx - boundary_size))) &
                           / (nakx - 2 * boundary_size)
               g0x(1, 1:boundary_size) = 0.0
               g0x(1, (nakx - boundary_size + 1):) = 0.0
            else
               g0x(1, :) = sum(g0x(1, :)) / nakx
            end if

            call transform_x2kx_unpadded(g0x, g0k)
            g1k = g1k + g0k
         end do

         phi_proj_stage(:, 1, it) = phi_proj_stage(:, 1, it) - g1k(1, :)
      end do
      deallocate (g0k, g1k, g0x)

   end subroutine add_adiabatic_response_radial

   !> Add radial variation of the Jacobian and gyroaveraing in the velocity integration of
   !> <g>, needed for radially global simulations
   subroutine add_radial_correction_int_species(g_in)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, is_idx
      use arrays_gyro_averages, only: aj0x, aj1x
      use geometry, only: dBdrho, bmag
      use arrays_store_useful, only: kperp2, dkperp2dr
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: vperp2
      use parameters_kxky_grid, only: nakx, naky
      use calculations_kxky, only: multiply_by_rho
      use parameters_multibox, only: ky_solve_radial
      use grids_species, only: spec

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(inout) :: g_in

      integer :: ivmu, iz, it, ia, imu, is, iky
      complex, dimension(:, :), allocatable :: g0k

      if (ky_solve_radial <= 0) return

      allocate (g0k(naky, nakx))

      ia = 1

      ! loop over super-index ivmu, which include vpa, mu and spec
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! is = species index
         is = is_idx(vmu_lo, ivmu)
         ! imu = mu index
         imu = imu_idx(vmu_lo, ivmu)

         ! loop over flux tubes in flux tube train
         do it = 1, ntubes
            ! loop over zed location within flux tube
            do iz = -nzgrid, nzgrid
               g0k = 0.0
               do iky = 1, min(ky_solve_radial, naky)
                  g0k(iky, :) = g_in(iky, :, iz, it, ivmu) &
                                * (-0.5 * aj1x(iky, :, iz, ivmu) / aj0x(iky, :, iz, ivmu) * (spec(is)%smz)**2 &
                                   * (kperp2(iky, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                   * (dkperp2dr(iky, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                   + dBdrho(iz) / bmag(ia, iz))

               end do
               !g0k(1,1) = 0.
               call multiply_by_rho(g0k)
               g_in(:, :, iz, it, ivmu) = g_in(:, :, iz, it, ivmu) + g0k
            end do
         end do
      end do

      deallocate (g0k)

   end subroutine add_radial_correction_int_species

   !> the following routine gets the correction in phi both from gyroaveraging and quasineutrality
   subroutine get_radial_correction(g, phi0, dist)

      use mp, only: proc0, mp_abort, sum_allreduce
      use stella_layouts, only: vmu_lo
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      use arrays_gyro_averages, only: aj0x, aj1x
      use parameters_physics, only: fphi
      use parameters_multibox, only: ky_solve_radial
      use geometry, only: dl_over_b, d_dl_over_b_drho, bmag, dBdrho
      use stella_layouts, only: imu_idx, is_idx
      use grids_z, only: nzgrid, ntubes
      use grids_velocity, only: integrate_species, vperp2
      use parameters_kxky_grid, only: nakx, naky
      use grids_kxky, only: rho_d_clamped
      use grids_kxky, only: zonal_mode
      use calculations_kxky, only: multiply_by_rho
      use grids_species, only: spec, has_electron_species
      use arrays_store_fields, only: phi_corr_QN, phi_corr_GA
      use arrays_store_fields, only: gamtot, dgamtotdr
      use arrays_store_fields, only: gamtot3, efac, efacp
      use arrays_store_useful, only: kperp2, dkperp2dr
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi0
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      character(*), intent(in) :: dist

      integer :: ikx, iky, ivmu, iz, it, ia, is, imu
      complex :: tmp
      real, dimension(:, :, :, :), allocatable :: gamtot_t
      complex, dimension(:, :, :, :), allocatable :: phi1
      complex, dimension(:, :, :), allocatable :: gyro_g
      complex, dimension(:, :), allocatable :: g0k, g1k, g1x

      ia = 1

      if (fphi > epsilon(0.0)) then
         allocate (gyro_g(naky, nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (g0k(naky, nakx))
         allocate (phi1(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi1 = 0.
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  is = is_idx(vmu_lo, ivmu)
                  imu = imu_idx(vmu_lo, ivmu)

                  g0k = g(:, :, iz, it, ivmu) &
                        * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) &
                           * (spec(is)%smz)**2 &
                           * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                           + dBdrho(iz) / bmag(ia, iz))

                  call gyro_average(g0k, iz, ivmu, gyro_g(:, :, ivmu))
               end do
               call integrate_species(gyro_g, iz, spec%z * spec%dens_psi0, phi1(:, :, iz, it), reduce_in=.false.)
            end do
         end do
         call sum_allreduce(phi1)

         !apply radial operator Xhat
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0k = phi1(:, :, iz, it) - dgamtotdr(:, :, iz) * phi0(:, :, iz, it)
               call multiply_by_rho(g0k)
               phi1(:, :, iz, it) = g0k
            end do
         end do

         if (dist == 'gbar') then
            allocate (gamtot_t(naky, nakx, -nzgrid:nzgrid, ntubes))
            gamtot_t = spread(gamtot, 4, ntubes)
            where (gamtot_t < epsilon(0.0))
               phi1 = 0.0
            elsewhere
               phi1 = phi1 / gamtot_t
            end where
            deallocate (gamtot_t)
         else if (dist == 'h') then
            if (proc0) write (*, *) 'dist option "h" not implemented in radial_correction. aborting'
            call mp_abort('dist option "h" in radial_correction. aborting')
         else
            if (proc0) write (*, *) 'unknown dist option in radial_correction. aborting'
            call mp_abort('unknown dist option in radial_correction. aborting')
         end if

         if (.not. has_electron_species(spec) .and. &
             adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
               if (dist == 'gbar') then
                  allocate (g1k(1, nakx))
                  allocate (g1x(1, nakx))
                  do it = 1, ntubes
                     do ikx = 1, nakx
                        g1k(1, ikx) = sum(phi0(1, ikx, :, it) &
                                          * (efacp * dl_over_b(ia, :) + efac * d_dl_over_b_drho(ia, :)))
                     end do
                     call transform_kx2x_unpadded(g1k, g1x)
                     g1x(1, :) = rho_d_clamped * g1x(1, :)
                     call transform_x2kx_unpadded(g1x, g1k)

                     do ikx = 1, nakx
                        phi1(1, ikx, :, it) = phi1(1, ikx, :, it) + g1k(1, ikx) / gamtot(1, ikx, :)
                        tmp = sum(dl_over_b(ia, :) * phi1(1, ikx, :, it))
                        phi1(1, ikx, :, it) = phi1(1, ikx, :, it) + gamtot3(ikx, :) * tmp
                     end do
                  end do
                  deallocate (g1k, g1x)
               else
                  if (proc0) write (*, *) 'unknown dist option in radial_correction. aborting'
                  call mp_abort('unknown dist option in radial_correction. aborting')
               end if
            end if
         end if

         !> collect quasineutrality corrections in wavenumber space
         phi_corr_QN = phi1

         !> zero out the ones we have already solved for using the full method
         do iky = 1, min(ky_solve_radial, naky)
            phi_corr_QN(iky, :, :, :) = 0.0
         end do

         deallocate (phi1)

         !> collect gyroaveraging corrections in wavenumber space
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call gyro_average_j1(phi0(:, :, iz, it), iz, ivmu, g0k)
                  g0k = -g0k * (spec(is)%smz)**2 &
                        * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                        * 0.5 * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))

                  call multiply_by_rho(g0k)
                  phi_corr_GA(:, :, iz, it, ivmu) = g0k
               end do
            end do
         end do

         deallocate (g0k)
         deallocate (gyro_g)

      end if

   end subroutine get_radial_correction

!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !============= INITALISE THE FIELDS FOR RADIALLY GLOBAL STELLA ==============
   !============================================================================
   subroutine init_fields_radial_variation
      use mp, only: proc0, job
      use mp, only: sum_allreduce
#ifdef ISO_C_BINDING
      use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_intptr_t
      use arrays_store_fields, only: qn_window, phi_shared
      use mp, only: sgproc0, curr_focus, sharedsubprocs
      use mp, only: scope, real_size, nbytes_real
      use mp, only: split_n_tasks, create_shared_memory_window
      use mpi
#endif
      use parameters_multibox, only: ky_solve_radial, ky_solve_real
      use grids_species, only: spec, has_electron_species, ion_species
      use calculations_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use grids_z, only: nzgrid, ntubes, nztot
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: akx
      use grids_kxky, only: zonal_mode, rho_d_clamped
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      use linear_solve, only: lu_decomposition, lu_inverse
      use multibox, only: init_mb_get_phi
      use arrays_store_fields, only: gamtot, dgamtotdr
      use arrays_store_fields, only: phi_solve, c_mat, theta
      use file_utils, only: runtype_option_switch, runtype_multibox

      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use arrays_gyro_averages, only: aj0v, aj1v
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use grids_velocity, only: vpa, vperp2, mu, nmu, nvpa
      use arrays_store_useful, only: kperp2, dkperp2dr
      use geometry, only: dBdrho, bmag
      use parameters_physics, only: tite, nine, beta

      use arrays_store_fields, only: efac, efacp
      use grids_velocity, only: integrate_vmu

      implicit none

      integer :: ikxkyz, it, is
      integer :: iz, ikx, iky, ia, zmi, naky_r
      real :: dum
      logical :: has_elec, adia_elec
#ifdef ISO_C_BINDING
      integer :: prior_focus, ierr
      integer :: counter, c_lo, c_hi
      integer(c_intptr_t):: cur_pos
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      complex, dimension(:), pointer :: phi_shared_temp
      type(c_ptr) :: cptr
#endif

      real, dimension(:, :), allocatable :: g0
      real, dimension(:), allocatable :: g1
      complex, dimension(:, :), allocatable :: g0k, g0x
      real :: wgt, tmp

      zm = 0
      ia = 1

      debug = debug .and. proc0

      call allocate_arrays_radial_variation

      allocate (g1(nmu))
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         it = it_idx(kxkyz_lo, ikxkyz)
         ! gamtot does not depend on flux tube index,
         ! so only compute for one flux tube index
         if (it /= 1) cycle
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         g1 = aj0v(:, ikxkyz) * aj1v(:, ikxkyz) * (spec(is)%smz)**2 &
              * (kperp2(iky, ikx, ia, iz) * vperp2(ia, iz, :) / bmag(ia, iz)**2) &
              * (dkperp2dr(iky, ikx, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
              / (1.0 - aj0v(:, ikxkyz)**2 + 100.*epsilon(0.0))

         g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa) &
              * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is) &
              * (-spec(is)%tprim * (spread(vpa**2, 2, nmu) + spread(vperp2(ia, iz, :), 1, nvpa) - 2.5) &
                 - spec(is)%fprim + (dBdrho(iz) / bmag(ia, iz)) * (1.0 - 2.0 * spread(mu, 1, nvpa) * bmag(ia, iz)) &
                 + spread(g1, 1, nvpa))
         wgt = spec(is)%z * spec(is)%z * spec(is)%dens / spec(is)%temp
         call integrate_vmu(g0, iz, tmp)
         dgamtotdr(iky, ikx, iz) = dgamtotdr(iky, ikx, iz) + tmp * wgt
      end do
      call sum_allreduce(dgamtotdr)
      deallocate (g1)

      if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
         dgamtotdr(1, 1, :) = 0.0
         zm = 1
      end if

      if (.not. has_electron_species(spec)) then
         efac = tite / nine * (spec(ion_species)%dens / spec(ion_species)%temp)
         efacp = efac * (spec(ion_species)%tprim - spec(ion_species)%fprim)
         dgamtotdr = dgamtotdr + efacp
      end if

      naky_r = min(naky, ky_solve_radial)

      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec .and. zonal_mode(1) &
                  .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      if (runtype_option_switch == runtype_multibox .and. job == 1 .and. ky_solve_real) then
         call init_mb_get_phi(has_elec, adia_elec, efac, efacp)
      elseif (runtype_option_switch /= runtype_multibox .or. (job == 1 .and. .not. ky_solve_real)) then
         allocate (g0k(1, nakx))
         allocate (g0x(1, nakx))

         if (.not. allocated(phi_solve)) allocate (phi_solve(naky_r, -nzgrid:nzgrid))
#ifdef ISO_C_BINDING
         prior_focus = curr_focus
         call scope(sharedsubprocs)
         !the following is to parallelize the calculation of QN for radial variation sims
         if (debug) write (*, *) 'fields::init_fields::phi_shared_init'
         if (phi_shared_window == MPI_WIN_NULL) then
            win_size = 0
            if (sgproc0) then
               win_size = int(naky * nakx * nztot * ntubes, MPI_ADDRESS_KIND) * 2 * real_size !complex size
            end if

            call create_shared_memory_window(win_size, phi_shared_window, cur_pos)

            cptr = transfer(cur_pos, cptr)

            if (.not. associated(phi_shared)) then
               ! associate array with lower bounds of 1
               call c_f_pointer(cptr, phi_shared_temp, (/naky * nakx * nztot * ntubes/))
               ! now get the correct bounds
               phi_shared(1:naky, 1:nakx, -nzgrid:nzgrid, 1:ntubes) => phi_shared_temp
            end if
            call mpi_win_fence(0, phi_shared_window, ierr)
         end if

         if (debug) write (*, *) 'fields::init_fields::qn_window_init'
         if (qn_window == MPI_WIN_NULL) then
            win_size = 0
            if (sgproc0) then
               win_size = int(nakx * nztot * naky_r, MPI_ADDRESS_KIND) * 4_MPI_ADDRESS_KIND &
                          + int(nakx**2 * nztot * naky_r, MPI_ADDRESS_KIND) * 2 * real_size !complex size
            end if

            call create_shared_memory_window(win_size, qn_window, cur_pos)

            !allocate the memory
            do iky = 1, naky_r
               zmi = 0
               if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
               do iz = -nzgrid, nzgrid
                  if (.not. associated(phi_solve(iky, iz)%zloc)) then
                     allocate (phi_solve(iky, iz)%zloc(nakx - zmi, nakx - zmi))
                     cptr = transfer(cur_pos, cptr)
                     call c_f_pointer(cptr, phi_solve(iky, iz)%zloc, (/nakx - zmi, nakx - zmi/))
                  end if
                  cur_pos = cur_pos + (nakx - zmi)**2 * 2 * nbytes_real
                  if (.not. associated(phi_solve(iky, iz)%idx)) then
                     cptr = transfer(cur_pos, cptr)
                     call c_f_pointer(cptr, phi_solve(iky, iz)%idx, (/nakx - zmi/))
                  end if
                  cur_pos = cur_pos + (nakx - zmi) * 4
               end do
            end do

            call mpi_win_fence(0, qn_window, ierr)
         end if

         call split_n_tasks(nztot * naky_r, c_lo, c_hi)

         call scope(prior_focus)
         counter = 0
#else
         do iky = 1, naky_r
            zmi = 0
            if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
            do iz = -nzgrid, nzgrid
               if (.not. associated(phi_solve(iky, iz)%zloc)) &
                  allocate (phi_solve(iky, iz)%zloc(nakx - zmi, nakx - zmi))
               if (.not. associated(phi_solve(iky, iz)%idx)) &
                  allocate (phi_solve(iky, iz)%idx(nakx - zmi))
            end do
         end do
#endif

         do iky = 1, naky_r
            zmi = 0
            if (iky == 1) zmi = zm !zero mode may or may not be included in matrix
            do iz = -nzgrid, nzgrid
#ifdef ISO_C_BINDING
               counter = counter + 1
               if ((counter >= c_lo) .and. (counter <= c_hi)) then
#endif
                  phi_solve(iky, iz)%zloc = 0.0
                  phi_solve(iky, iz)%idx = 0
                  do ikx = 1 + zmi, nakx
                     g0k(1, :) = 0.0
                     g0k(1, ikx) = dgamtotdr(iky, ikx, iz)

                     call transform_kx2x_unpadded(g0k, g0x)
                     g0x(1, :) = rho_d_clamped * g0x(1, :)
                     call transform_x2kx_unpadded(g0x, g0k)

                     !row column
                     phi_solve(iky, iz)%zloc(:, ikx - zmi) = g0k(1, (1 + zmi):)
                     phi_solve(iky, iz)%zloc(ikx - zmi, ikx - zmi) = phi_solve(iky, iz)%zloc(ikx - zmi, ikx - zmi) &
                                                                     + gamtot(iky, ikx, iz)
                  end do

                  call lu_decomposition(phi_solve(iky, iz)%zloc, phi_solve(iky, iz)%idx, dum)
#ifdef ISO_C_BINDING
               end if
#endif
            end do
         end do

         if (adia_elec) then
            if (.not. allocated(c_mat)) allocate (c_mat(nakx, nakx)); 
            if (.not. allocated(theta)) allocate (theta(nakx, nakx, -nzgrid:nzgrid)); 
            !get C
            do ikx = 1, nakx
               g0k(1, :) = 0.0
               g0k(1, ikx) = 1.0

               call transform_kx2x_unpadded(g0k, g0x)
               g0x(1, :) = (efac + efacp * rho_d_clamped) * g0x(1, :)
               call transform_x2kx_unpadded(g0x, g0k)

               !row column
               c_mat(:, ikx) = g0k(1, :)
            end do

            !get Theta
            do iz = -nzgrid, nzgrid

               !get Theta
               do ikx = 1, nakx
                  g0k(1, :) = 0.0
                  g0k(1, ikx) = dgamtotdr(1, ikx, iz) - efacp

                  call transform_kx2x_unpadded(g0k, g0x)
                  g0x(1, :) = rho_d_clamped * g0x(1, :)
                  call transform_x2kx_unpadded(g0x, g0k)

                  !row column
                  theta(:, ikx, iz) = g0k(1, :)
                  theta(ikx, ikx, iz) = theta(ikx, ikx, iz) + gamtot(1, ikx, iz) - efac
               end do
            end do
         end if
         deallocate (g0k, g0x)
      end if

   end subroutine init_fields_radial_variation

   !============================================================================
   !======================= ALLOCATE ARRAYS FOR EM FIELDS ======================
   !============================================================================

   subroutine allocate_arrays_radial_variation

      use parameters_physics, only: radial_variation
      use arrays_store_fields, only: phi_corr_QN, phi_corr_GA
      use arrays_store_fields, only: apar_corr_QN, apar_corr_GA
      use arrays_store_fields, only: dgamtotdr
      use grids_z, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo
      use parameters_kxky_grid, only: naky, nakx

      implicit none

      if (.not. allocated(phi_corr_QN) .and. radial_variation) then
         allocate (phi_corr_QN(naky, nakx, -nzgrid:nzgrid, ntubes))
         phi_corr_QN = 0.
      end if
      if (.not. allocated(phi_corr_GA) .and. radial_variation) then
         allocate (phi_corr_GA(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         phi_corr_GA = 0.
      end if
      if (.not. allocated(apar_corr_QN) .and. radial_variation) then
         !allocate (apar_corr(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (apar_corr_QN(1, 1, 1, 1))
         apar_corr_QN = 0.
      end if
      if (.not. allocated(apar_corr_GA) .and. radial_variation) then
         !allocate (apar_corr(naky,nakx,-nzgrid:nzgrid,ntubes,vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (apar_corr_GA(1, 1, 1, 1, 1))
         apar_corr_GA = 0.
      end if

      if (.not. allocated(dgamtotdr)) allocate (dgamtotdr(naky, nakx, -nzgrid:nzgrid)); dgamtotdr = 0.
      if (.not. allocated(dgamtotdr)) allocate (dgamtotdr(1, 1, 1)); dgamtotdr = 0.

   end subroutine allocate_arrays_radial_variation

   !============================================================================
   !==================== FINISH THE RADIALLY GLOBAL FIELDS =====================
   !============================================================================
   subroutine finish_radial_fields

      use arrays_store_fields, only: phi_corr_QN, phi_corr_GA
      use arrays_store_fields, only: apar_corr_QN, apar_corr_GA
      use arrays_store_fields, only: dgamtotdr
      use arrays_store_fields, only: c_mat, theta
#ifdef ISO_C_BINDING
      use arrays_store_fields, only: qn_window
      use mpi
#endif

      implicit none

#ifdef ISO_C_BINDING
      integer ierr
#endif

#ifdef ISO_C_BINDING
      if (phi_shared_window /= MPI_WIN_NULL) call mpi_win_free(phi_shared_window, ierr)
      if (qn_window /= MPI_WIN_NULL) then
         call mpi_win_free(qn_window, ierr)
      end if
#endif

      if (allocated(c_mat)) deallocate (c_mat)
      if (allocated(theta)) deallocate (theta)

      if (allocated(phi_corr_QN)) deallocate (phi_corr_QN)
      if (allocated(phi_corr_GA)) deallocate (phi_corr_GA)

      if (allocated(apar_corr_QN)) deallocate (apar_corr_QN)
      if (allocated(apar_corr_GA)) deallocate (apar_corr_GA)

      if (allocated(dgamtotdr)) deallocate (dgamtotdr)

   end subroutine finish_radial_fields

end module fields_radial_variation
