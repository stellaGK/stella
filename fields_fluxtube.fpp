module fields_fluxtube

   implicit none

   public :: init_fields_fluxtube
   public :: advance_fields_fluxtube
   public :: get_fields, get_fields_vmulo
   public :: get_fields_by_spec, get_fields_by_spec_idx

   !> NB: probably want to put these in a folder so EM and radial stella can use/modify them

   private

   integer :: zm
   logical :: debug = .false.

   !> init_fields_fluxtube allocates and fills arrays needed during main time advance
   !> loop for the field solve for flux tube simulations

contains

   !> ######################################################################################
   !> Initialise fields arrays for EM stella
   !> ######################################################################################

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

      ! could move these array allocations to allocate_arrays to clean up code
      if (.not. allocated(gamtot3)) then
         if (.not. has_electron_species(spec) &
             .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            allocate (gamtot3(nakx, -nzgrid:nzgrid)); gamtot3 = 0.
         else
            allocate (gamtot3(1, 1)); gamtot3 = 0.
         end if
      end if

      if (fphi > epsilon(0.0)) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! gamtot does not depend on flux tube index,
            ! so only compute for one flux tube index
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa)
            if (.not. maxwellian_normalization) then
               g0 = g0 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            end if
            wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            gamtot(iky, ikx, iz) = gamtot(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(gamtot)

         gamtot_h = sum(spec%z * spec%z * spec%dens / spec%temp)

         ! avoid divide by zero when kx=ky=0
         ! do not evolve this mode, so value is irrelevant
         if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
            gamtot(1, 1, :) = 0.0
            zm = 1
         end if

         if (.not. has_electron_species(spec)) then
            efac = tite / nine * (spec(ion_species)%dens / spec(ion_species)%temp)
            efacp = efac * (spec(ion_species)%tprim - spec(ion_species)%fprim)
            gamtot = gamtot + efac
            gamtot_h = gamtot_h + efac
            if (adiabatic_option_switch == adiabatic_option_fieldlineavg) then
               if (zonal_mode(1)) then
                  gamtot3_h = efac / (sum(spec%zt * spec%z * spec%dens))
                  do ikx = 1, nakx
                     ! avoid divide by zero for kx=ky=0 mode,
                     ! which we do not need anyway
                     !if (abs(akx(ikx)) < epsilon(0.)) cycle
                     tmp = 1./efac - sum(dl_over_b(ia, :) / gamtot(1, ikx, :))
                     gamtot3(ikx, :) = 1./(gamtot(1, ikx, :) * tmp)
                  end do
                  if (akx(1) < epsilon(0.)) then
                     gamtot3(1, :) = 0.0
                  end if
               end if
            end if
         end if

         deallocate (g0)

      end if

   end subroutine init_fields_fluxtube

   !> ######################################################################################
   !> This is to advance fields for fluxtube stella using Quasineutrality in the main code
   !> ######################################################################################
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
         !> first gather (vpa,mu) onto processor for v-space operations
         !> v-space operations are field solve, dg/dvpa, and collisions
         if (debug) write (*, *) 'dist_fn::advance_stella::scatter'
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         call scatter(kxkyz2vmu, g, gvmu)
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         !> given gvmu with vpa and mu local, calculate the corresponding fields
         if (debug) write (*, *) 'dist_fn::advance_stella::get_fields'
         call get_fields(gvmu, phi, apar, dist)
      else
         call get_fields_vmulo(g, phi, apar, dist)
      end if
   end subroutine advance_fields_fluxtube

   !> TODO-GA: remove apar from this and make it only needed for EM stella
   subroutine get_fields(g, phi, apar, dist, skip_fsa)

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

      if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_kxkyzlo'

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
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            wgt = spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
         end do
         deallocate (g0)
         call sum_allreduce(phi)
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         call get_phi(phi, dist, skip_fsa_local)

      end if

   end subroutine get_fields

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

      if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_vmulo'

      phi = 0.
      if (fphi > epsilon(0.0)) then
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         ! gyroaverage the distribution function g at each phase space location
         call gyro_average(g, g_gyro)

         ! <g> requires modification if radial profile variation is included
         if (radial_variation) call add_radial_correction_int_species(g_gyro)

         ! integrate <g> over velocity space and sum over species
         !> store result in phi, which will be further modified below to account for polarization term
         if (debug) write (*, *) 'dist_fn::advance_stella::sum_all_reduce'
         call integrate_species(g_gyro, spec%z * spec%dens_psi0, phi)

         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')

         call get_phi(phi, dist, skip_fsa_local)

      end if

   end subroutine get_fields_vmulo

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

         ! divide <g> by sum_s (\Gamma_0s-1) Zs^2*e*ns/Ts to get phi
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

   !> ######################################################################################
   !> Routines Used for Collision Operator
   !> ######################################################################################
   subroutine get_fields_by_spec(g, fld, skip_fsa)

      use mp, only: sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use gyro_averages, only: gyro_average
      use run_parameters, only: fphi
      use stella_geometry, only: dl_over_b
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: nvpa, nmu
      use vpamu_grids, only: integrate_vmu
      use kt_grids, only: nakx
      use kt_grids, only: zonal_mode
      use species, only: spec, nspec, has_electron_species
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg

      use fields_arrays, only: gamtot3_h, gamtot_h

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
      logical, optional, intent(in) :: skip_fsa

      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local
      complex, dimension(nspec) :: tmp

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      if (debug) write (*, *) 'dist_fn::advance_stella::get_fields_by_spec'

      ia = 1

      fld = 0.
      if (fphi > epsilon(0.0)) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            wgt = spec(is)%z * spec(is)%dens_psi0
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
            g0 = g0 * wgt
            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
         end do
         call sum_allreduce(fld)

         fld = fld / gamtot_h

         if (.not. has_electron_species(spec) .and. (.not. skip_fsa_local) .and. &
             adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
               do ikx = 1, nakx
                  do it = 1, ntubes
                     do is = 1, nspec
                        tmp(is) = sum(dl_over_b(ia, :) * fld(1, ikx, :, it, is))
                        fld(1, ikx, :, it, is) = fld(1, ikx, :, it, is) + tmp(is) * gamtot3_h
                     end do
                  end do
               end do
            end if
         end if

         deallocate (g0)
      end if

   end subroutine get_fields_by_spec

   subroutine get_fields_by_spec_idx(isa, g, fld)

      ! apply phi_isa[ ] to all species indices contained in g
      ! ie get phi_isa[g_is1], phi_isa[g_is2], phi_isa[g_is3] ...

      use mp, only: sum_allreduce
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      use gyro_averages, only: gyro_average
      use run_parameters, only: fphi
      use stella_geometry, only: dl_over_b, bmag
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: vperp2, nvpa, nmu
      use vpamu_grids, only: integrate_vmu
      use kt_grids, only: nakx
      use kt_grids, only: zonal_mode
      use species, only: spec, nspec, has_electron_species
      use physics_flags, only: adiabatic_option_switch
      use physics_flags, only: adiabatic_option_fieldlineavg
      use dist_fn_arrays, only: kperp2
      use spfunc, only: j0

      use fields_arrays, only: gamtot_h, gamtot3_h
      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld
      integer, intent(in) :: isa

      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia, imu
      complex, dimension(nspec) :: tmp
      real :: wgt
      real :: arg

      ia = 1

      fld = 0.
      if (fphi > epsilon(0.0)) then
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            wgt = spec(isa)%z * spec(isa)%dens
            do imu = 1, nmu
               ! AVB: changed this for use of j0, check
               arg = spec(isa)%bess_fac * spec(isa)%smz_psi0 * sqrt(vperp2(ia, iz, imu) * kperp2(iky, ikx, ia, iz)) / bmag(ia, iz)
               g0(:, imu) = g(:, imu, ikxkyz) * j0(arg) ! AVB: gyroaverage
            end do
            g0 = g0 * wgt
            call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
         end do
         call sum_allreduce(fld)

         fld = fld / gamtot_h

         if (.not. has_electron_species(spec) .and. &
             adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            if (zonal_mode(1)) then
               do ikx = 1, nakx
                  do it = 1, ntubes
                     do is = 1, nspec
                        tmp(is) = sum(dl_over_b(ia, :) * fld(1, ikx, :, it, is))
                        fld(1, ikx, :, it, is) = fld(1, ikx, :, it, is) + tmp(is) * gamtot3_h
                     end do
                  end do
               end do
            end if
         end if

         deallocate (g0)
      end if

   end subroutine get_fields_by_spec_idx

end module fields_fluxtube
