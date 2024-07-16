module fields_fluxtube

   implicit none

   !> Advance Routine
   public :: advance_fields_fluxtube
   !> Initialise Routine
   public :: init_fields_fluxtube
   public :: get_fields_kxkyzlo

   private

   integer :: zm
   logical :: debug = .false.

contains

!###############################################################################
!############################## ADVANCE FIELDS #################################
!###############################################################################

   !============================================================================
   !============================== ADVANCE FIELDS ==============================
   !============================================================================
   !> This calls the appropriate routines needed to advance phi in the main code
   !> when using fluxtube stella, depending on the distribution (i.e. if the
   !> information is parallelised over (kx,ky,z) or (vpa,mu) ).
   !> Note that Apar and Bpar are only advanced when using EM so these are in
   !> fields_electromagnetic.fpp
   !============================================================================
   subroutine advance_fields_fluxtube(g, phi, apar, dist)

      use mp, only: proc0
      use stella_layouts, only: vmu_lo
      use job_manage, only: time_message
      use redistribute, only: scatter
      use dist_fn_arrays, only: gvmu
      use zgrid, only: nzgrid
      use dist_redistribute, only: kxkyz2vmu
      use run_parameters, only: fields_kxkyz
      use physics_flags, only: full_flux_surface
      use fields_arrays, only: time_field_solve

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      character(*), intent(in) :: dist

      !> fields_kxkyz = F is the default
      if (fields_kxkyz) then
         !> First gather (vpa,mu) onto processor for v-space operations
         !> v-space operations are field solve, dg/dvpa, and collisions
         if (debug) write (*, *) 'fields::fields_fluxtube::advance_fields_fluxtube::scatter'
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         call scatter(kxkyz2vmu, g, gvmu)
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         !> Given gvmu with vpa and mu local, calculate the corresponding fields
         if (debug) write (*, *) 'fields::fields_fluxtube::advance_fields_fluxtube::get_fields'
         call get_fields_kxkyzlo(gvmu, phi, apar, dist)
      else
         call get_fields_vmulo(g, phi, apar, dist)
      end if
   end subroutine advance_fields_fluxtube

   !============================================================================
   !=========================== ADVANCE FIELDS VMULO ===========================
   !============================================================================
   !> If we are parallelising over (vpa,mu) then this subroutine is called
   !> This is the more common version used compared with parallelising over (kx,ky,z)
   !> Here we calculate:
   !>    sum_s int dv (J0 * g)
   !> and then call get_phi which divides this by the appropriate gamtot factor
   !============================================================================
   !> TODO-GA: remove apar from this and make it only needed for EM stella
   subroutine get_fields_vmulo(g, phi, apar, dist, skip_fsa)

      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use gyro_averages, only: gyro_average
      use run_parameters, only: fphi
      use physics_flags, only: radial_variation
      use dist_fn_arrays, only: g_gyro
      use zgrid, only: nzgrid
      use vpamu_grids, only: integrate_species
      use species, only: spec

      use fields_radial_variation, only: add_radial_correction_int_species
      use fields_arrays, only: time_field_solve

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist

      logical :: skip_fsa_local

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      if (debug) write (*, *) 'fields::fields_fluxtube::get_fields_vmulo'

      phi = 0.
      if (fphi > epsilon(0.0)) then
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         !> Gyroaverage the distribution function g at each phase space location
         call gyro_average(g, g_gyro)

         !> TODO-GA: Move this somehow??
         ! If radial profile variation is included then <g> requires modification
         if (radial_variation) call add_radial_correction_int_species(g_gyro)

         !> Integrate <g> over velocity space and sum over species
         !> Store result in phi, which will be further modified below to account for polarization term
         if (debug) write (*, *) 'fields::fields_fluxtube::get_fields_vmulo::integrate_species'
         call integrate_species(g_gyro, spec%z * spec%dens_psi0, phi)

         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         !> Divide this by the appropriate factor that appears in Quasineutrality
         call get_phi(phi, dist, skip_fsa_local)

      end if

   end subroutine get_fields_vmulo

   !============================================================================
   !=========================== ADVANCE FIELDS KXKYZ ===========================
   !============================================================================
   !> If we are parallelising over (kx,ky,z) then this subroutine is called
   !> This is the less common version used compared with parallelising over (vpa, mu)
   !> Here we calculate:
   !>    sum_s int dv (J0 * g)
   !> and then call get_phi which divides this by the appropriate gamtot factor
   !============================================================================
   !> TODO-GA: remove apar from this and make it only needed for EM stella
   subroutine get_fields_kxkyzlo(g, phi, apar, dist, skip_fsa)

      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use gyro_averages, only: gyro_average
      use run_parameters, only: fphi
      use physics_parameters, only: beta
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: vpa
      use vpamu_grids, only: integrate_vmu
      use species, only: spec
      use fields_arrays, only: time_field_solve

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist
      complex :: tmp

      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      if (debug) write (*, *) 'fields::fields_fluxtube::get_fields_kxkyzlo'

      ia = 1

      phi = 0.
      if (fphi > epsilon(0.0)) then
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            !> Gyroaverage the distribution function g at each phase space location
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)

            !> Integrate <g> over velocity space
            call integrate_vmu(g0, iz, tmp)

            !> Sum over species
            !> Store result in phi, which will be further modified below to account for polarization term
            wgt = spec(is)%z * spec(is)%dens_psi0
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
         end do
         deallocate (g0)

         !> Gather information that is spread over different processors onto one and then
         !> broadcast to all processors
         call sum_allreduce(phi)
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         !> Divide this by the appropriate factor that appears in Quasineutrality
         call get_phi(phi, dist, skip_fsa_local)

      end if

   end subroutine get_fields_kxkyzlo

   !============================================================================
   !================================ UPDATE PHI ================================
   !============================================================================
   !> The 'phi' variable passed in is
   !>    sum_s int dv (J0 * g)
   !> This routine divides by the appropriate gamtot factor depending on if we
   !> have kinetic or adiabatic electrons, and also on whether we are using 'g'
   !> or 'h' as our distribution function that we are evolving
   !============================================================================
   subroutine get_phi(phi, dist, skip_fsa)

      use mp, only: proc0, mp_abort, job
      use job_manage, only: time_message
      use physics_flags, only: radial_variation
      use run_parameters, only: ky_solve_radial, ky_solve_real
      use zgrid, only: nzgrid, ntubes
      use stella_geometry, only: dl_over_b
      use kt_grids, only: nakx, naky, zonal_mode
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use species, only: spec, has_electron_species
      use multibox, only: mb_get_phi
      use fields_arrays, only: gamtot, gamtot3, gamtot_h, gamtot3_h
      use file_utils, only: runtype_option_switch, runtype_multibox
      use fields_arrays, only: time_field_solve

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
      if (dist == 'h') then
         phi = phi / gamtot_h
      else if (dist == 'gbar') then

         !> Divide <g> by sum_s (1 - Gamma_0) Z^2*n/T to get phi
         allocate (gamtot_t(naky, nakx, -nzgrid:nzgrid, ntubes))
         gamtot_t = spread(gamtot, 4, ntubes)
         where (gamtot_t < epsilon(0.0))
            phi = 0.0
         elsewhere
            phi = phi / gamtot_t
         end where
         deallocate (gamtot_t)
      else
         if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
         call mp_abort('unknown dist option in get_fields. aborting')
         return
      end if

      if (any(gamtot(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')

      !> Now handle adiabatic electrons if needed
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
         else if (dist == 'gbar') then
            do ikx = 1, nakx
               do it = 1, ntubes
                  tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                  phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * gamtot3(ikx, :)
               end do
            end do
         else
            if (proc0) write (*, *) 'unknown dist option in get_fields. aborting'
            call mp_abort('unknown dist option in get_fields. aborting')
         end if
      end if
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'get_phi_adia_elec')

   end subroutine get_phi

!###############################################################################
!################################## INITALISE ##################################
!###############################################################################

   !============================================================================
   !============================= INITALISE ARRAYS =============================
   !============================================================================
   !> Initialise arrays that are needed during the main time advance loop in the
   !> field solve for the flux tube simulations only
   !> These are initialised once and then used throughout the rest of the simulation
   !>    gamtot = 1 - gamma0 = Z^2*n/T * int e^(-v^2) * (1 - J0^2) dv
   !> If using adiabatic electrons then this factor is modified and we use gamtot3
   !> which includes the Boltmann response
   !============================================================================
   subroutine init_fields_fluxtube

      use mp, only: proc0, sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use gyro_averages, only: aj0v
      use run_parameters, only: fphi
      use run_parameters, only: maxwellian_normalization
      use physics_parameters, only: tite, nine, beta
      use species, only: spec, has_electron_species, ion_species
      use stella_geometry, only: dl_over_b, bmag
      use zgrid, only: nzgrid
      use vpamu_grids, only: nvpa, nmu, mu
      use vpamu_grids, only: vpa, vperp2
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: integrate_vmu
      use kt_grids, only: naky, nakx, akx
      use kt_grids, only: zonal_mode
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use fields_arrays, only: gamtot, gamtot3
      use fields_arrays, only: gamtot_h, gamtot3_h, efac, efacp

      implicit none

      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      real :: tmp, wgt
      real, dimension(:, :), allocatable :: g0

      debug = debug .and. proc0

      ia = 1
      zm = 0

      !> Initialise gamtot
      if (fphi > epsilon(0.0)) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            !> Gamtot does not depend on flux tube index, so only compute for one flux tube index
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            !> Calculate (1- j0^2) for each kx,ky,z
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa)
            !> TODO-GA: remove maxwellian_normalisation flag
            if (.not. maxwellian_normalization) then
               g0 = g0 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            end if

            !> Integrate over vpa mu with weight = Z^2*n/T
            wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            gamtot(iky, ikx, iz) = gamtot(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(gamtot)

         gamtot_h = sum(spec%z * spec%z * spec%dens / spec%temp)

         !> Avoid divide by zero when kx=ky=0
         !> We do not evolve this mode, so the value is irrelevant
         if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
            gamtot(1, 1, :) = 0.0
            zm = 1
         end if

         !> Initialise gamtot3
         if (.not. has_electron_species(spec)) then
            efac = tite / nine * (spec(ion_species)%dens / spec(ion_species)%temp)
            efacp = efac * (spec(ion_species)%tprim - spec(ion_species)%fprim)
            gamtot = gamtot + efac
            gamtot_h = gamtot_h + efac
            if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
               if (zonal_mode(1)) then
                  gamtot3_h = efac / (sum(spec%zt * spec%z * spec%dens))
                  do ikx = 1, nakx
                     tmp = 1./efac - sum(dl_over_b(ia, :) / gamtot(1, ikx, :))
                     gamtot3(ikx, :) = 1./(gamtot(1, ikx, :) * tmp)
                  end do
                  !> Avoid divide by zero for kx=ky=0 mode, which we do not need anyway
                  if (akx(1) < epsilon(0.)) then
                     gamtot3(1, :) = 0.0
                  end if
               end if
            end if
         end if

         deallocate (g0)

      end if

   end subroutine init_fields_fluxtube

end module fields_fluxtube
