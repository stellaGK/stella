! Module for advancing and initialising the fields for a fluxtube simulation
module fields_fluxtube

   use debug_flags, only: debug => fields_fluxtube_debug
   implicit none

   ! Advance fields for fluxtube routines
   public :: advance_fields_fluxtube
   public :: get_fields_fluxtube
   ! Initialise Routine
   public :: init_fields_fluxtube   

   private

   interface get_fields_fluxtube
      module procedure get_fields_fluxtube_vmlo
      module procedure get_fields_fluxtube_kxkyzlo
   end interface

   integer :: zm

contains

!###############################################################################
!############################## ADVANCE FIELDS #################################
!###############################################################################

   !============================================================================
   !============================== ADVANCE FIELDS ==============================
   !============================================================================
   ! This calls the appropriate routines needed to advance phi in the main code
   ! when using fluxtube stella, depending on the distribution (i.e. if the
   ! information is parallelised over (kx,ky,z) or (vpa,mu) ).
   ! Note that Apar and Bpar are only advanced when using EM so these are in
   ! fields_electromagnetic.fpp
   !============================================================================
   subroutine advance_fields_fluxtube (g, phi, apar, bpar, dist, skip_fsa)

      use mp, only: proc0
      use job_manage, only: time_message
      ! Layouts
      use stella_layouts, only: vmu_lo
      use redistribute, only: scatter
      use calculations_redistribute, only: kxkyz2vmu
      ! Arrays
      use arrays_store_distribution_fn, only: gvmu
      use arrays_store_useful, only: time_field_solve
      ! Parameters
      use stella_layouts, only: fields_kxkyz
      ! Grids
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: phi, apar, bpar
      character(*), intent(in) :: dist

      logical, optional, intent(in) :: skip_fsa
      !-------------------------------------------------------------------------
      ! Note that fields_kxkyz = F is the default
      if (.not. fields_kxkyz) then 
         if (debug) write (*, *) 'fields_fluxtube::advance_fields_fluxtube::get_fields_fluxtube_vmlo'
         ! This will call get_fields_fluxtube_vmlo
         call get_fields_fluxtube(g, phi, apar, bpar, dist, skip_fsa)
      else if (fields_kxkyz) then 
         ! Note that to use this option it has to be specified by the user
         ! First gather (vpa,mu) onto processor for v-space operations
         ! v-space operations are field solve, dg/dvpa, and collisions
         if (debug) write (*, *) 'fields_fluxtube::advance_fields_fluxtube::scatter'
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         call scatter(kxkyz2vmu, g, gvmu)
         if (proc0) call time_message(.false., time_field_solve(:, 2), ' fields_redist')
         ! Given gvmu with vpa and mu local, calculate the corresponding fields
         ! This will call get_fields_fluxtube_kxkyzlo
         if (debug) write (*, *) 'fields_fluxtube::advance_fields_fluxtube::get_fields_fluxtube_kxkyzlo'
         call get_fields_fluxtube(gvmu, phi, apar, bpar, dist, skip_fsa)
      end if 

   end subroutine advance_fields_fluxtube

   !============================================================================
   !=========================== ADVANCE FIELDS VMULO ===========================
   !============================================================================
   ! If we are parallelising over (vpa,mu) then this subroutine is called
   ! This is the more common version used compared with parallelising over 
   ! (kx,ky,z) and is the default for stella.
   ! Here we calculate:
   !    sum_s int dv (J0 * g)
   ! and then call get_phi which divides this by the appropriate gamtot factor
   !============================================================================
   ! TODO-GA: remove apar from this and make it only needed for EM stella
   subroutine get_fields_fluxtube_vmlo(g, phi, apar, bpar, dist, skip_fsa)

      use mp, only: mp_abort, proc0
      use job_manage, only: time_message
      ! Layouts
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx
      ! Arrays
      use arrays_store_distribution_fn, only: g_scratch
      use arrays_store_useful, only: time_field_solve
      ! Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: radial_variation
      use parameters_physics, only: fphi
      ! Grids
      use grids_z, only: nzgrid
      use grids_species, only: spec
      use calculations_velocity_integrals, only: integrate_species
      ! Calculations
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      ! Routines from other field modules
      use fields_electromagnetic, only: get_fields_electromagnetic
      use fields_radial_variation, only: add_radial_correction_int_species, get_phi_for_radial
      
      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: phi, apar, bpar
      character(*), intent(in) :: dist 

      logical, optional, intent(in) :: skip_fsa
      logical :: skip_fsa_local
      !-------------------------------------------------------------------------
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa
      if (debug) write (*, *) 'fields_fluxtube::get_fields_fluxtube_vmulo'

      phi = 0.
      ! Note that this advances phi for the electrostatic, fluxtube case.
      ! If electromagnetic effects are included then phi will be advanced below
      if (fphi > epsilon(0.0) .and. .not. include_bpar) then 
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         ! First gyroaverage the distribution function g at each phase space location
         ! and store this as g_scratch = <g> = J_0 g in k-space
         call gyro_average(g, g_scratch)
         ! If we are allowing for Radial Variation then we must modify <g>.
         ! This is done in the fields_radialvariation module, but is not needed
         ! for standard stella (as these are false by default).
         if (radial_variation) call add_radial_correction_int_species(g_scratch)
         ! Next, integrate <g> over velocity space and sum over species and store
         ! the result in phi, which will be further modified below to account for polarization term
         if (debug) write (*, *) 'fields_fluxtube::get_fields_fluxtube_vmulo::integrate_species_phi'
         call integrate_species(g_scratch, spec%z * spec%dens_psi0, phi)
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         ! Get phi routine inverts the appropriate 'Gamma_0' factor that multiplies the <phi> terms
         ! onto the RHS of the equation. The exact form of this factor depends on the simulation
         ! e.g. we may have adiabatic electrons, kinetic electrons etc. and also if Radial Varation 
         ! effects are included (these are false by default)
         if (.not. radial_variation) then 
            ! This is the inversion routine for electrostatic, fluxtube stella
            if (debug) write (*, *) 'fields_fluxtube::get_fields_fluxtube_vmulo::get_phi'
            call get_phi(phi, dist, skip_fsa_local)
         else if (radial_variation) then 
            ! If we are allowing for radial variation this factor changes further
            if (debug) write (*, *) 'fields_fluxtube::get_fields_fluxtube_vmulo::get_phi_for_radial'
            call get_phi_for_radial (phi, dist , skip_fsa)
         end if 
      end if 

      ! If we have electromagnetic effects we need to add the terms due to A_parallel and B_parallel 
      ! and also need to advance phi including electromagnetic effects
      if (include_apar .or. include_bpar) then
         apar = 0.
         bpar = 0.
         if (debug) write (*, *) 'fields_fluxtube::get_fields_fluxtube_vmulo::get_fields_electromagnetic'
         call get_fields_electromagnetic(g, phi, apar, bpar, dist)
      end if 

   end subroutine get_fields_fluxtube_vmlo

   !============================================================================
   !=========================== ADVANCE FIELDS KXKYZ ===========================
   !============================================================================
   ! If we are parallelising over (kx,ky,z) then this subroutine is called
   ! This is the less common version used compared with parallelising over 
   ! (vpa, mu). This is NOT the default routine for stella, and the flag 
   ! <fields_kxkyz = .True.> must be set to use this routine. 
   ! Here we calculate:
   !    sum_s int dv (J0 * g)
   ! and then call get_phi which divides this by the appropriate gamtot factor
   !============================================================================
   ! TODO-GA: remove apar from this and make it only needed for EM stella
   subroutine get_fields_fluxtube_kxkyzlo(g, phi, apar, bpar, dist, skip_fsa)

      use mp, only: proc0
      use mp, only: sum_allreduce, mp_abort
      use job_manage, only: time_message
      ! Layouts
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      ! Arrays
      use arrays_store_useful, only: time_field_solve
      ! Parameters
      use parameters_physics, only: include_apar, include_bpar
      use parameters_physics, only: fphi
      use parameters_physics, only: radial_variation
      ! Grids
      use grids_velocity, only: nvpa, nmu
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_species, only: spec
      use grids_z, only: nzgrid
      ! Calculations
      use calculations_gyro_averages, only: gyro_average, gyro_average_j1
      ! Routines from other field modules
      use fields_electromagnetic, only: get_fields_electromagnetic
      use fields_radial_variation, only: get_phi_for_radial

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(inout) :: phi, apar, bpar
      logical, optional, intent(in) :: skip_fsa
      character(*), intent(in) :: dist
      complex :: tmp

      real :: wgt
      complex, dimension(:, :), allocatable :: g0
      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      logical :: skip_fsa_local
      !-------------------------------------------------------------------------
      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      if (debug) write (*, *) 'fields_fluxtube::get_fields_fluxtube_kxkyzlo'

      ia = 1
      phi = 0.

      ! Note that this advances phi for the electrostatic, fluxtube case.
      ! If electromagnetic effects are included then phi will be advanced below
      if (fphi > epsilon(0.0) .and. .not. include_bpar) then
         if (proc0) call time_message(.false., time_field_solve(:, 3), ' int_dv_g')
         if (debug) write (*, *) 'fields::get_fields_fluxtube_kxkyzlo_int_dv_g'
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iky = iky_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            ! First gyroaverage the distribution function g at each phase space location
            ! and store this as g0 = <g> = J_0 g in k-space
            call gyro_average(g(:, :, ikxkyz), ikxkyz, g0)
             ! Next, integrate <g> over velocity space and sum over species and store
            ! the result in phi, which will be further modified below to account for polarization term
            wgt = spec(is)%z * spec(is)%dens_psi0
            call integrate_vmu(g0, iz, tmp)
            phi(iky, ikx, iz, it) = phi(iky, ikx, iz, it) + wgt * tmp
         end do
         deallocate (g0)
         call sum_allreduce(phi)

         if (.not. radial_variation) then 
            if (debug) write (*, *) 'fields_fluxtube::get_fields_fluxtube_kxkyzlo::get_phi'
            call get_phi(phi, dist, skip_fsa_local)
         else if (radial_variation) then
            ! If we are allowing for radial variation this factor changes further
            if (debug) write (*, *) 'fields_fluxtube::get_fields_fluxtube_kxkyzlo::get_phi_for_radial'
            call get_phi_for_radial (phi, dist , skip_fsa)
         end if
      end if

      ! If we have electromagnetic effects we need to add the terms due to A_parallel and B_parallel 
      ! and also need to advance phi including electromagnetic effects
      if (include_apar .or. include_bpar) then 
         bpar = 0.
         apar = 0.
         if (debug) write (*, *) 'fields_fluxtube::advance_fields_kxkyzlo::get_fields_electromagnetic'
         call get_fields_electromagnetic(g, phi, apar, bpar, dist)
      end if 

   end subroutine get_fields_fluxtube_kxkyzlo

   !============================================================================
   !================================ UPDATE PHI ================================
   !============================================================================
   ! The 'phi' variable passed in is:
   !    sum_s int dv (J0 * g)
   ! This routine divides by the appropriate gamtot factor depending on if we
   ! have kinetic or adiabatic electrons, and also on whether we are using 'g'
   ! or 'h' as our distribution function that we are evolving.
   ! Note that this routine is only called in the Electrostatic, Fluxtube case. 
   !============================================================================
   subroutine get_phi(phi, dist, skip_fsa)

      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use multibox, only: mb_get_phi
      ! Arrays
      use arrays_store_useful, only: gamtot, gamtot3, gamtot_h, gamtot3_h
      use arrays_store_useful, only: time_field_solve
      ! Parameters
      use grids_kxky, only: nakx, naky
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: zonal_mode
      use grids_species, only: spec, has_electron_species
      ! Geometry 
      use geometry, only: dl_over_b

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi
      logical, optional, intent(in) :: skip_fsa
      real, dimension(:, :, :, :), allocatable :: gamtot_t
      integer :: ia, it, ikx
      complex :: tmp
      logical :: skip_fsa_local
      logical :: has_elec, adia_elec

      character(*), intent(in) :: dist
      !-------------------------------------------------------------------------
      if (debug) write (*, *) 'fields_fluxtube::get_phi'

      skip_fsa_local = .false.
      if (present(skip_fsa)) skip_fsa_local = skip_fsa

      ia = 1
      has_elec = has_electron_species(spec)
      adia_elec = .not. has_elec &
                  .and. adiabatic_option_switch == adiabatic_option_fieldlineavg

      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')
      if (dist == 'h') then
         if (debug) write(*, *) 'fields_fluxtube::get_phi::dist==h'
         phi = phi / gamtot_h
      else if (dist == 'g' .or. dist == 'gbar') then
         if (debug) write(*, *) 'fields_fluxtube::get_phi::dist==gbar'
         ! Divide <g> by sum_s (1 - Gamma_0) Z^2*n/T to get phi
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

      ! The kx = ky = 0.0 mode is not evolved by stella so make sure this term
      ! is set to zero.
      if (debug) write(*, *) 'fields_fluxtube::get_phi::set kxky=0.0 to zero'
      if (any(gamtot(1, 1, :) < epsilon(0.))) phi(1, 1, :, :) = 0.0
      if (proc0) call time_message(.false., time_field_solve(:, 4), ' get_phi')

      ! Now handle adiabatic electrons only if needed.
      if (proc0) call time_message(.false., time_field_solve(:, 5), 'get_phi_adia_elec')
      if (adia_elec .and. zonal_mode(1) .and. .not. skip_fsa_local) then
         if (debug) write (*, *) 'dist_fn::advance_stella::adiabatic_electrons'
         if (dist == 'h') then
            if (debug) write(*, *) 'fields_fluxtube::get_phi::dist==h adiabatic'
            do it = 1, ntubes
               do ikx = 1, nakx
                  tmp = sum(dl_over_b(ia, :) * phi(1, ikx, :, it))
                  phi(1, ikx, :, it) = phi(1, ikx, :, it) + tmp * gamtot3_h
               end do
            end do
         else if (dist == 'g' .or. dist == 'gbar') then
            if (debug) write(*, *) 'fields_fluxtube::get_phi::dist==gbar adaiabtic'
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
   ! Initialise arrays that are needed during the main time advance loop in the
   ! field solve for the flux tube simulations only
   ! These are initialised once and then used throughout the rest of the simulation
   !    gamtot = 1 - gamma0 = Z^2*n/T * int e^(-v^2) * (1 - J0^2) dv
   ! If using adiabatic electrons then this factor is modified and we use gamtot3
   ! which includes the Boltmann response
   !============================================================================
   subroutine init_fields_fluxtube

      use mp, only: proc0, sum_allreduce
      ! Layouts
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, onlY: iz_idx, it_idx, ikx_idx, iky_idx, is_idx
      ! Arrays
      use arrays_store_useful, only: gamtot, gamtot3
      use arrays_store_useful, only: gamtot_h, gamtot3_h, efac, efacp
      ! Parameters
      use parameters_physics, only: fphi
      use parameters_numerical, only: maxwellian_normalization
      use parameters_physics, only: tite, nine
      use grids_kxky, only: nakx
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg
      ! Grids
      use grids_velocity, only: nvpa, nmu
      use grids_velocity, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use calculations_velocity_integrals, only: integrate_vmu
      use grids_species, only: spec, has_electron_species, ion_species
      use grids_kxky, only: zonal_mode, akx
      ! Calculations
      use arrays_gyro_averages, only: aj0v
      ! Geometry
      use geometry, only: dl_over_b

      implicit none

      integer :: ikxkyz, iz, it, ikx, iky, is, ia
      real :: tmp, wgt
      real, dimension(:, :), allocatable :: g0
      !-------------------------------------------------------------------------
      debug = debug .and. proc0

      ia = 1
      zm = 0

      ! Initialise gamtot
      if (fphi > epsilon(0.0)) then
         if (debug) write(*, *) 'fields_fluxtube::init_fields_fluxtube::init_gamtot'
         allocate (g0(nvpa, nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            it = it_idx(kxkyz_lo, ikxkyz)
            ! Gamtot does not depend on flux tube index, so only compute for one flux tube index
            if (it /= 1) cycle
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            ! Calculate (1- j0^2) for each kx,ky,z
            g0 = spread((1.0 - aj0v(:, ikxkyz)**2), 1, nvpa)
            ! TODO-GA: remove maxwellian_normalisation flag
            if (.not. maxwellian_normalization) then
               g0 = g0 * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * maxwell_fac(is)
            end if

            ! Integrate over vpa mu with weight = Z^2*n/T
            wgt = spec(is)%z * spec(is)%z * spec(is)%dens_psi0 / spec(is)%temp
            call integrate_vmu(g0, iz, tmp)
            gamtot(iky, ikx, iz) = gamtot(iky, ikx, iz) + tmp * wgt
         end do
         call sum_allreduce(gamtot)

         gamtot_h = sum(spec%z * spec%z * spec%dens / spec%temp)

         ! Avoid divide by zero when kx=ky=0
         ! We do not evolve this mode, so the value is irrelevant
         if (zonal_mode(1) .and. akx(1) < epsilon(0.) .and. has_electron_species(spec)) then
            gamtot(1, 1, :) = 0.0
            zm = 1
         end if

         ! Initialise gamtot3
         if (.not. has_electron_species(spec)) then
            if (debug) write(*, *) 'fields_fluxtube::init_fields_fluxtube::init_gamtot3'
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
                  ! Avoid divide by zero for kx=ky=0 mode, which we do not need anyway
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
