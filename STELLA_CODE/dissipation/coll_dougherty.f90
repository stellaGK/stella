module coll_dougherty

   implicit none

   public :: read_parameters_dougherty
   public :: init_collisions_dougherty
   public :: finish_collisions_dougherty
   public :: advance_collisions_dougherty_explicit
   public :: advance_collisions_dougherty_implicit

   private

   logical :: vpa_operator, mu_operator
   logical :: momentum_conservation, energy_conservation

   integer :: nresponse_vpa = 1
   integer :: nresponse_mu = 1

   real, dimension(:, :), allocatable :: aa_vpa, bb_vpa, cc_vpa
   real, dimension(:, :, :), allocatable :: aa_mu, cc_mu
   real, dimension(:, :), allocatable :: bb_mu
   complex, dimension(:, :, :), allocatable :: vpadiff_response
   integer, dimension(:, :), allocatable :: vpadiff_idx
   complex, dimension(:, :, :), allocatable :: mudiff_response
   integer, dimension(:, :), allocatable :: mudiff_idx

   complex, dimension(:, :, :), allocatable :: vpadiff_zf_response
   integer, dimension(:, :), allocatable :: vpadiff_zf_idx
   complex, dimension(:, :, :), allocatable :: mudiff_zf_response
   integer, dimension(:, :), allocatable :: mudiff_zf_idx

   logical :: dougherty_initialized = .false.

contains

   subroutine read_parameters_dougherty

      use mp, only: broadcast
      use namelist_dissipation, only: read_namelist_collisions_dougherty

      implicit none

      call read_namelist_collisions_dougherty(momentum_conservation, energy_conservation, vpa_operator, mu_operator) 

      call broadcast(momentum_conservation)
      call broadcast(energy_conservation)
      call broadcast(vpa_operator)
      call broadcast(mu_operator)

   end subroutine read_parameters_dougherty

   subroutine init_collisions_dougherty(collisions_implicit, cfl_dt_vpadiff, cfl_dt_mudiff)

      use species, only: spec, nspec
      use velocity_grids, only: dvpa, dmu, mu, nmu
      use geometry, only: bmag

      implicit none

      logical, intent(in) :: collisions_implicit
      real, intent(out) :: cfl_dt_vpadiff, cfl_dt_mudiff

      integer :: is
      real :: vnew_max

      if (dougherty_initialized) return
      dougherty_initialized = .true.

      if (collisions_implicit) then
         if (vpa_operator) then
            call init_vpadiff_matrix
            call init_vpadiff_conserve
         end if
         if (mu_operator) then
            call init_mudiff_matrix
            call init_mudiff_conserve
         end if
      else
         vnew_max = 0.0
         do is = 1, nspec
            vnew_max = max(vnew_max, maxval(spec(is)%vnew))
         end do
         cfl_dt_vpadiff = 2.0 * dvpa**2 / vnew_max
         cfl_dt_mudiff = minval(bmag) / (vnew_max * maxval(mu(2:) / dmu(:nmu - 1)**2))
      end if
   end subroutine init_collisions_dougherty

   subroutine init_vpadiff_matrix

      use stella_time, only: code_dt
      use species, only: nspec, spec
      use velocity_grids, only: dvpa, vpa, nvpa
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
      use geometry, only: bmag
      use arrays_dist_fn, only: kperp2

      implicit none

      integer :: ikxkyz, iky, ikx, iz, is
      integer :: ia

      if (.not. allocated(aa_vpa)) allocate (aa_vpa(nvpa, nspec))
!    if (.not.allocated(bb_vpa)) allocate (bb_vpa(nvpa,nspec))
      if (.not. allocated(bb_vpa)) allocate (bb_vpa(nvpa, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      if (.not. allocated(cc_vpa)) allocate (cc_vpa(nvpa, nspec))

      ! deal with boundary points (BC is f(vpa)=0 beyond +/- vpa_max)
      aa_vpa(1, :) = 0.0; cc_vpa(nvpa, :) = 0.0
      ! 2nd order centered differences for d/dvpa (1/2 dh/dvpa + vpa h)
      do is = 1, nspec
         aa_vpa(2:, is) = -code_dt * spec(is)%vnew(is) * 0.5 * (1.0 / dvpa - vpa(:nvpa - 1)) / dvpa
!       bb_vpa(:,is) = 1.0+code_dt*spec(is)%vnew(is)/dvpa**2
         cc_vpa(:nvpa - 1, is) = -code_dt * spec(is)%vnew(is) * 0.5 * (1.0 / dvpa + vpa(2:)) / dvpa
      end do

      ia = 1
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         bb_vpa(:, ikxkyz) = 1.0 + code_dt * spec(is)%vnew(is) &
                             * (0.25 * kperp2(iky, ikx, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2 + 1./dvpa**2)
      end do

   end subroutine init_vpadiff_matrix

   subroutine init_mudiff_matrix

      use stella_time, only: code_dt
      use species, only: nspec, spec
      use z_grid, only: nzgrid
      use geometry, only: bmag
      use velocity_grids, only: dmu, nmu
      use velocity_grids, only: dmu_cell, mu_cell, wgts_mu_bare
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, is_idx
      use arrays_dist_fn, only: kperp2

      implicit none

      integer :: ikxkyz, iky, ikx, iz, is
      integer :: ia
      ! TMP FOR TESTING -- MAB
!    integer :: imu

      if (.not. allocated(aa_mu)) allocate (aa_mu(-nzgrid:nzgrid, nmu, nspec))
      if (.not. allocated(bb_mu)) allocate (bb_mu(nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      if (.not. allocated(cc_mu)) allocate (cc_mu(-nzgrid:nzgrid, nmu, nspec))

      ia = 1

      ! deal with boundary points (BC is f(mu)=0 beyond mu_max and collision operator vanishes for mu -> 0)
      aa_mu(:, 1, :) = 0.0; cc_mu(:, nmu, :) = 0.0
      ! 2nd order centered differences for dt * nu * d/dmu (mu/B*dh/dmu + 2*mu*h)
      do is = 1, nspec
         do iz = -nzgrid, nzgrid
            aa_mu(iz, 2:, is) = -code_dt * spec(is)%vnew(is) * mu_cell(:nmu - 1) * (1.0 / (bmag(ia, iz) * dmu) - 1.0) / wgts_mu_bare(2:)
            cc_mu(iz, :nmu - 1, is) = -code_dt * spec(is)%vnew(is) * mu_cell(:nmu - 1) * (1.0 / (bmag(ia, iz) * dmu) + 1.0) / wgts_mu_bare(:nmu - 1)
         end do
      end do

      !1st derivative here is  2 d/dmu( mu h) -> (mu(i+1/2)h(i+1/2) - mu(i-1/2)h(i-1/2))/(mu(i+1/2)-mu(i-1/2))
      !where h(i+1/2) = 0.5*[h(i+1)+h(i)] and mu(i+1/2) = 0.5*(mu(i+1)+mu(i))
      !2nd derivative here is d/dmu ( mu dh/dmu ) = (mu(i+1/2)h'(i+1/2) - mu(i-1/2)h'(i-1/2))/(mu(i+1/2)-mu(i-1/2))
      !where h'(i+1/2) = (h(i+1)-h(i))/(mu(i+1)-mu(i))
      !left  endpoint i=1   -> mu(i-1/2) = 0
      !right endpoint i=nmu -> h(i+1/2) = h'(i+1/2) = 0 (these are the two boundary conditions)
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         bb_mu(1, ikxkyz) = 1.0 + code_dt * spec(is)%vnew(is) &
                            * (0.25 * kperp2(iky, ikx, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2 &
                               + dmu_cell(1) * (1.0 / (dmu(1) * bmag(ia, iz)) - 1.0) / wgts_mu_bare(1))
         bb_mu(2:nmu - 1, ikxkyz) = 1.0 + code_dt * spec(is)%vnew(is) &
                                    * (0.25 * kperp2(iky, ikx, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2 &
                                       + (mu_cell(2:nmu - 1) / dmu(2:) + mu_cell(:nmu - 2) / dmu(:nmu - 2)) &
                                       / (wgts_mu_bare(2:nmu - 1) * bmag(ia, iz)) - dmu_cell(2:nmu - 1) / wgts_mu_bare(2:nmu - 1))
         bb_mu(nmu, ikxkyz) = 1.0 + code_dt * spec(is)%vnew(is) &
                              * (0.25 * kperp2(iky, ikx, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2 &
                                 + mu_cell(nmu - 1) * (1.0 / (dmu(nmu - 1) * bmag(ia, iz)) + 1.0) / wgts_mu_bare(nmu))
      end do

   end subroutine init_mudiff_matrix

   subroutine init_vpadiff_conserve

      use mp, only: sum_allreduce
      use finite_differences, only: tridag
      use linear_solve, only: lu_decomposition, lu_inverse
      use stella_time, only: code_dt
      use species, only: nspec, spec, has_electron_species
      use z_grid, only: nzgrid, ntubes
      use velocity_grids, only: ztmax, maxwell_vpa, maxwell_mu
      use velocity_grids, only: nmu, vpa, vperp2
      use velocity_grids, only: set_vpa_weights
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: zonal_mode
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use geometry, only: dl_over_b
      use arrays_dist_fn, only: gvmu
      use gyro_averages, only: aj0v
      use fields_fluxtube, only: get_fields_fluxtube
      use fields_collisions, only: get_fields_by_spec
      use arrays_fields, only: efac, gamtot_h
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg

      implicit none

      integer :: ikxkyz, iky, ikx, iz, it, is, ia
      integer :: imu
      integer :: idx
      logical :: conservative_wgts
      real :: dum2
      complex, dimension(:, :, :, :), allocatable :: dum1
      complex, dimension(:, :, :, :), allocatable :: dum3
      complex, dimension(:, :, :, :, :), allocatable :: field
      complex, dimension(:, :), allocatable :: temp_mat

      ia = 1

      nresponse_vpa = 1
      if (momentum_conservation) nresponse_vpa = nresponse_vpa + nspec
      if (energy_conservation) nresponse_vpa = nresponse_vpa + nspec

      if (.not. allocated(vpadiff_response)) then
         allocate (vpadiff_response(nresponse_vpa, nresponse_vpa, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         vpadiff_response = 0.
         allocate (vpadiff_idx(nresponse_vpa, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      end if

      if (.not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         if (.not. allocated(vpadiff_zf_response)) then
            allocate (vpadiff_zf_response(nresponse_vpa, nresponse_vpa, nakx))
            vpadiff_zf_response = 0.
            allocate (vpadiff_zf_idx(nresponse_vpa, nakx))
         end if
      end if

      allocate (dum1(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dum3(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      ! set wgts to be equally spaced to ensure exact conservation properties
      conservative_wgts = .true.
      call set_vpa_weights(conservative_wgts)

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            gvmu(:, imu, ikxkyz) = ztmax(:, is) * maxwell_mu(1, iz, imu, is) * aj0v(imu, ikxkyz)
            call tridag(1, aa_vpa(:, is), bb_vpa(:, ikxkyz), cc_vpa(:, is), gvmu(:, imu, ikxkyz))
         end do
      end do

      ! gvmu contains dhs/dphi
      ! for phi equation, need 1-P[dhs/dphi]
      ! for upar equations, need -Us[dhs/dphi]
      ! for energy conservation, need -Qs[dhs/dphi]
      call get_fields_fluxtube(gvmu, field(:, :, :, :, 1), dum1, dum3, dist='h', skip_fsa=.true.)

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         vpadiff_response(1, 1, ikxkyz) = 1.0 - field(iky, ikx, iz, it, 1)
      end do
      idx = 2
      if (momentum_conservation) then
         call get_upar(gvmu, field)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            vpadiff_response(idx:idx + nspec - 1, 1, ikxkyz) = -field(iky, ikx, iz, it, :)
         end do
         idx = idx + nspec
      end if
      if (energy_conservation) then
         call get_temp(gvmu, field)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            vpadiff_response(idx:idx + nspec - 1, 1, ikxkyz) = -field(iky, ikx, iz, it, :)
         end do
      end if
      idx = 2

      if (momentum_conservation) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            do imu = 1, nmu
               gvmu(:, imu, ikxkyz) = 2.*code_dt * spec(is)%vnew(is) * vpa * aj0v(imu, ikxkyz) * maxwell_vpa(:, is) * maxwell_mu(1, iz, imu, is)
               call tridag(1, aa_vpa(:, is), bb_vpa(:, ikxkyz), cc_vpa(:, is), gvmu(:, imu, ikxkyz))
            end do
         end do
         ! gvmu now contains dhs/dupars
         ! need to get -Ps[dhs/dupars] for phi equation
         ! need to get 1-Us[dhs/dupars] for momentum conservation
         ! need to get -Qs[dhs/dupars] for energy conservation
         call get_fields_by_spec(gvmu, field, skip_fsa=.true.)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            vpadiff_response(1, idx:idx + nspec - 1, ikxkyz) = -field(iky, ikx, iz, it, :)
         end do

         call get_upar(gvmu, field)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            do is = 1, nspec
               vpadiff_response(idx + is - 1, idx + is - 1, ikxkyz) = 1.0 - field(iky, ikx, iz, it, is)
            end do
         end do

         if (energy_conservation) then
            call get_temp(gvmu, field)
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iky = iky_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iz = iz_idx(kxkyz_lo, ikxkyz)
               it = it_idx(kxkyz_lo, ikxkyz)
               do is = 1, nspec
                  vpadiff_response(idx + is + nspec - 1, idx + is - 1, ikxkyz) = -field(iky, ikx, iz, it, is)
               end do
            end do
         end if
         idx = idx + nspec
      end if

      if (energy_conservation) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            do imu = 1, nmu
               gvmu(:, imu, ikxkyz) = 2.*code_dt * spec(is)%vnew(is) * (vpa**2 + vperp2(1, iz, imu) - 1.5) &
                                      * aj0v(imu, ikxkyz) * maxwell_vpa(:, is) * maxwell_mu(1, iz, imu, is)
               call tridag(1, aa_vpa(:, is), bb_vpa(:, ikxkyz), cc_vpa(:, is), gvmu(:, imu, ikxkyz))
            end do
         end do
         ! gvmu now contains dhs/dQs
         ! need to get -Ps[dhs/dQs] for phi equation
         ! need to get 1-Us[dhs/dQs] for momentum conservation
         ! need to get -Qs[dhs/dQs] for energy conservation
         call get_fields_by_spec(gvmu, field, skip_fsa=.true.)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            vpadiff_response(1, idx:idx + nspec - 1, ikxkyz) = -field(iky, ikx, iz, it, :)
         end do

         if (momentum_conservation) then
            call get_upar(gvmu, field)
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iky = iky_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iz = iz_idx(kxkyz_lo, ikxkyz)
               it = it_idx(kxkyz_lo, ikxkyz)
               do is = 1, nspec
                  vpadiff_response(idx + is - 1 - nspec, idx + is - 1, ikxkyz) = -field(iky, ikx, iz, it, is)
               end do
            end do
         end if

         call get_temp(gvmu, field)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            do is = 1, nspec
               vpadiff_response(idx + is - 1, idx + is - 1, ikxkyz) = 1.0 - field(iky, ikx, iz, it, is)
            end do
         end do
      end if

      ! now get LU decomposition for vpadiff_response
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         call lu_decomposition(vpadiff_response(:, :, ikxkyz), vpadiff_idx(:, ikxkyz), dum2)
      end do

      ! if electrons are adiabatic, compute the matrices for the flux-surface-average
      if (.not. has_electron_species(spec) .and. zonal_mode(1) &
          .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         allocate (temp_mat(nresponse_vpa, nresponse_vpa))
         vpadiff_zf_response = 0.0
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            if (iky /= 1 .or. it /= 1) cycle

            !calculate inverse of vpadiff_response
            call lu_inverse(vpadiff_response(:, :, ikxkyz), vpadiff_idx(:, ikxkyz), temp_mat)

            !calculate -inv(vpadiff_response).Q, where Q has a single entry
            do idx = 1, nresponse_vpa
               vpadiff_zf_response(idx, 1, ikx) = vpadiff_zf_response(idx, 1, ikx) &
                                                  - temp_mat(idx, 1) * (efac / gamtot_h) * dl_over_b(ia, iz)
            end do
         end do

         !finish the flux surface average
         call sum_allreduce(vpadiff_zf_response)

         !calculate 1 - fsa(inv(vpadiff_response).Q)
         do idx = 1, nresponse_vpa
            vpadiff_zf_response(idx, idx, :) = vpadiff_zf_response(idx, idx, :) + 1.0
         end do

         do ikx = 1, nakx
            call lu_decomposition(vpadiff_zf_response(:, :, ikx), vpadiff_zf_idx(:, ikx), dum2)
         end do

         deallocate (temp_mat)
      end if

      ! reset wgts to default setting
      conservative_wgts = .false.
      call set_vpa_weights(conservative_wgts)

      deallocate (dum1, dum3, field)

   end subroutine init_vpadiff_conserve

   subroutine init_mudiff_conserve

      use mp, only: sum_allreduce
      use finite_differences, only: tridag
      use linear_solve, only: lu_decomposition, lu_inverse
      use stella_time, only: code_dt
      use species, only: nspec, spec, has_electron_species
      use z_grid, only: nzgrid, ntubes
      use velocity_grids, only: ztmax, maxwell_vpa, maxwell_mu
      use velocity_grids, only: nvpa, vpa, vperp2
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: zonal_mode
      use geometry, only: dl_over_b, bmag
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use arrays_dist_fn, only: gvmu, kperp2
      use gyro_averages, only: aj0v, aj1v
      use fields_fluxtube, only: get_fields_fluxtube
      use fields_collisions, only: get_fields_by_spec
      use arrays_fields, only: efac, gamtot_h
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg

      implicit none

      integer :: ikxkyz, iky, ikx, iz, it, is, ia
      integer :: iv
      integer :: idx
      real :: dum2
      complex, dimension(:, :), allocatable :: temp_mat
      complex, dimension(:, :, :, :), allocatable :: dum1, dum3
      complex, dimension(:, :, :, :, :), allocatable :: field

      nresponse_mu = 1
      if (momentum_conservation) nresponse_mu = nresponse_mu + nspec
      if (energy_conservation) nresponse_mu = nresponse_mu + nspec

      if (.not. allocated(mudiff_response)) then
         allocate (mudiff_response(nresponse_mu, nresponse_mu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         mudiff_response = 0.
         allocate (mudiff_idx(nresponse_mu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      end if

      if (.not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         if (.not. allocated(mudiff_zf_response)) then
            allocate (mudiff_zf_response(nresponse_mu, nresponse_mu, nakx))
            mudiff_zf_response = 0.
            allocate (mudiff_zf_idx(nresponse_mu, nakx))
         end if
      end if

      allocate (dum1(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dum3(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do iv = 1, nvpa
            gvmu(iv, :, ikxkyz) = ztmax(iv, is) * maxwell_mu(1, iz, :, is) * aj0v(:, ikxkyz)
            call tridag(1, aa_mu(iz, :, is), bb_mu(:, ikxkyz), cc_mu(iz, :, is), gvmu(iv, :, ikxkyz))
         end do
      end do

      ! gvmu contains dhs/dphi
      ! for phi equation, need 1-P[dhs/dphi]
      ! for uperp equations, need -Us[dhs/dphi]
      ! for energy conservation, need -Qs[dhs/dphi]
      call get_fields_fluxtube(gvmu, field(:, :, :, :, 1), dum1, dum3, dist='h', skip_fsa=.true.)

      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         mudiff_response(1, 1, ikxkyz) = 1.0 - field(iky, ikx, iz, it, 1)
      end do
      idx = 2
      if (momentum_conservation) then
         call get_uperp(gvmu, field)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            mudiff_response(idx:idx + nspec - 1, 1, ikxkyz) = -field(iky, ikx, iz, it, :)
         end do
         idx = idx + nspec
      end if
      if (energy_conservation) then
         call get_temp_mu(gvmu, field)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            mudiff_response(idx:idx + nspec - 1, 1, ikxkyz) = -field(iky, ikx, iz, it, :)
         end do
      end if
      idx = 2

      if (momentum_conservation) then
         ia = 1
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            do iv = 1, nvpa
               gvmu(iv, :, ikxkyz) = 2.*code_dt * spec(is)%vnew(is) * kperp2(iky, ikx, ia, iz) * vperp2(ia, iz, :) &
                                     * (spec(is)%smz / bmag(ia, iz))**2 * aj1v(:, ikxkyz) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, :, is)
               call tridag(1, aa_mu(iz, :, is), bb_mu(:, ikxkyz), cc_mu(iz, :, is), gvmu(iv, :, ikxkyz))
            end do
         end do
         ! gvmu now contains dhs/dupars
         ! need to get -Ps[dhs/dupars] for phi equation
         ! need to get 1-Us[dhs/dupars] for momentum conservation
         ! need to get -Qs[dhs/dupars] for energy conservation
         call get_fields_by_spec(gvmu, field, skip_fsa=.true.)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            mudiff_response(1, idx:idx + nspec - 1, ikxkyz) = -field(iky, ikx, iz, it, :)
         end do

         call get_uperp(gvmu, field)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            do is = 1, nspec
               mudiff_response(idx + is - 1, idx + is - 1, ikxkyz) = 1.0 - field(iky, ikx, iz, it, is)
            end do
         end do

         if (energy_conservation) then
            call get_temp_mu(gvmu, field)
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iky = iky_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iz = iz_idx(kxkyz_lo, ikxkyz)
               it = it_idx(kxkyz_lo, ikxkyz)
               do is = 1, nspec
                  mudiff_response(idx + is + nspec - 1, idx + is - 1, ikxkyz) = -field(iky, ikx, iz, it, is)
               end do
            end do
         end if
         idx = idx + nspec
      end if

      if (energy_conservation) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            do iv = 1, nvpa
               gvmu(iv, :, ikxkyz) = 2.0 * code_dt * spec(is)%vnew(is) * (vpa(iv)**2 + vperp2(1, iz, :) - 1.5) &
                                     * aj0v(:, ikxkyz) * maxwell_vpa(iv, is) * maxwell_mu(1, iz, :, is)
               call tridag(1, aa_mu(iz, :, is), bb_mu(:, ikxkyz), cc_mu(iz, :, is), gvmu(iv, :, ikxkyz))
            end do
         end do
         ! gvmu now contains dhs/dQs
         ! need to get -Ps[dhs/dQs] for phi equation
         ! need to get 1-Us[dhs/dQs] for momentum conservation
         ! need to get -Qs[dhs/dQs] for energy conservation
         call get_fields_by_spec(gvmu, field, skip_fsa=.true.)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            mudiff_response(1, idx:idx + nspec - 1, ikxkyz) = -field(iky, ikx, iz, it, :)
         end do

         if (momentum_conservation) then
            call get_uperp(gvmu, field)
            do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
               iky = iky_idx(kxkyz_lo, ikxkyz)
               ikx = ikx_idx(kxkyz_lo, ikxkyz)
               iz = iz_idx(kxkyz_lo, ikxkyz)
               it = it_idx(kxkyz_lo, ikxkyz)
               do is = 1, nspec
                  mudiff_response(idx + is - 1 - nspec, idx + is - 1, ikxkyz) = -field(iky, ikx, iz, it, is)
               end do
            end do
         end if

         call get_temp_mu(gvmu, field)
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            do is = 1, nspec
               mudiff_response(idx + is - 1, idx + is - 1, ikxkyz) = 1.0 - field(iky, ikx, iz, it, is)
            end do
         end do
      end if

      ! now get LU decomposition for mudiff_response
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         call lu_decomposition(mudiff_response(:, :, ikxkyz), mudiff_idx(:, ikxkyz), dum2)
      end do

      ! if electrons are adiabatic, compute the matrices for the flux-surface-average
      if (.not. has_electron_species(spec) .and. zonal_mode(1) &
          .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         allocate (temp_mat(nresponse_mu, nresponse_mu))
         mudiff_zf_response = 0.0
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            if (iky /= 1 .or. it /= 1) cycle

            !calculate inverse of mudiff_response
            call lu_inverse(mudiff_response(:, :, ikxkyz), mudiff_idx(:, ikxkyz), temp_mat)

            !calculate -inv(mudiff_response).Q, where Q has a single entry
            do idx = 1, nresponse_mu
               mudiff_zf_response(idx, 1, ikx) = mudiff_zf_response(idx, 1, ikx) &
                                                 - temp_mat(idx, 1) * (efac / gamtot_h) * dl_over_b(ia, iz)
            end do
         end do

         !finish the flux surface average
         call sum_allreduce(mudiff_zf_response)

         !calculate 1 - fsa(inv(mudiff_response).Q)
         do idx = 1, nresponse_mu
            mudiff_zf_response(idx, idx, :) = mudiff_zf_response(idx, idx, :) + 1.0
         end do

         do ikx = 1, nakx
            call lu_decomposition(mudiff_zf_response(:, :, ikx), mudiff_zf_idx(:, ikx), dum2)
         end do

         deallocate (temp_mat)
      end if

      deallocate (dum1, dum3, field)

   end subroutine init_mudiff_conserve

   subroutine finish_collisions_dougherty

      implicit none

      call finish_vpadiff_matrix
      call finish_mudiff_matrix
      call finish_vpadiff_response
      call finish_mudiff_response

      dougherty_initialized = .false.

   end subroutine finish_collisions_dougherty

   subroutine finish_vpadiff_matrix

      implicit none

      if (allocated(aa_vpa)) deallocate (aa_vpa)
      if (allocated(bb_vpa)) deallocate (bb_vpa)
      if (allocated(cc_vpa)) deallocate (cc_vpa)

   end subroutine finish_vpadiff_matrix

   subroutine finish_mudiff_matrix

      implicit none

      if (allocated(aa_mu)) deallocate (aa_mu)
      if (allocated(bb_mu)) deallocate (bb_mu)
      if (allocated(cc_mu)) deallocate (cc_mu)

   end subroutine finish_mudiff_matrix

   subroutine finish_vpadiff_response

      implicit none

      if (allocated(vpadiff_response)) deallocate (vpadiff_response)
      if (allocated(vpadiff_idx)) deallocate (vpadiff_idx)

   end subroutine finish_vpadiff_response

   subroutine finish_mudiff_response

      implicit none

      if (allocated(mudiff_response)) deallocate (mudiff_response)
      if (allocated(mudiff_idx)) deallocate (mudiff_idx)

   end subroutine finish_mudiff_response

   subroutine get_upar(g, fld)

      use mp, only: sum_allreduce
      use z_grid, only: nzgrid
      use velocity_grids, only: integrate_vmu
      use velocity_grids, only: nvpa, nmu
      use velocity_grids, only: vpa
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use gyro_averages, only: aj0v

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld

      integer :: ikxkyz, iky, ikx, iz, it, is
      complex, dimension(:, :), allocatable :: g0

      allocate (g0(nvpa, nmu))

      fld = 0.
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         g0 = g(:, :, ikxkyz) * spread(vpa, 2, nmu) * spread(aj0v(:, ikxkyz), 1, nvpa)
         call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
      end do
      deallocate (g0)

      call sum_allreduce(fld)

   end subroutine get_upar

   subroutine get_uperp(g, fld)

      use mp, only: sum_allreduce
      use z_grid, only: nzgrid
      use velocity_grids, only: integrate_vmu
      use velocity_grids, only: nvpa, nmu
      use velocity_grids, only: vperp2
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use gyro_averages, only: aj1v

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld

      integer :: ikxkyz, iky, ikx, iz, it, is
      complex, dimension(:, :), allocatable :: g0

      allocate (g0(nvpa, nmu))

      fld = 0.
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
!        g0 = 2.0*g(:,:,ikxkyz)*spread((vperp2(1,iz,:)-0.5)*aj1v(:,ikxkyz),1,nvpa)
         g0 = g(:, :, ikxkyz) * spread((vperp2(1, iz, :) - 0.5) * aj1v(:, ikxkyz), 1, nvpa)
         call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
      end do
      deallocate (g0)

      call sum_allreduce(fld)

   end subroutine get_uperp

   subroutine get_temp(g, fld)

      use mp, only: sum_allreduce
      use z_grid, only: nzgrid
      use velocity_grids, only: integrate_vmu
      use velocity_grids, only: nvpa, nmu
      use velocity_grids, only: vpa
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use gyro_averages, only: aj0v

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld

      integer :: ikxkyz, iky, ikx, iz, it, is
      complex, dimension(:, :), allocatable :: g0

      allocate (g0(nvpa, nmu))

      fld = 0.
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         g0 = g(:, :, ikxkyz) * (spread(vpa**2, 2, nmu) - 0.5) &
              * spread(aj0v(:, ikxkyz), 1, nvpa) / 1.5
         call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
      end do
      deallocate (g0)

      call sum_allreduce(fld)

   end subroutine get_temp

   subroutine get_temp_mu(g, fld)

      use mp, only: sum_allreduce
      use z_grid, only: nzgrid
      use velocity_grids, only: integrate_vmu
      use velocity_grids, only: nvpa, nmu, vperp2
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use gyro_averages, only: aj0v

      implicit none

      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: fld

      integer :: ikxkyz, iky, ikx, iz, it, is
      complex, dimension(:, :), allocatable :: g0

      allocate (g0(nvpa, nmu))

      fld = 0.
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         g0 = g(:, :, ikxkyz) * (spread(vperp2(1, iz, :), 1, nvpa) - 1.0) &
              * spread(aj0v(:, ikxkyz), 1, nvpa) / 1.5
         call integrate_vmu(g0, iz, fld(iky, ikx, iz, it, is))
      end do
      deallocate (g0)

      call sum_allreduce(fld)

   end subroutine get_temp_mu

   subroutine advance_collisions_dougherty_explicit(g, phi, bpar, gke_rhs, time_collisions)

      use mp, only: proc0, mp_abort
      use job_manage, only: time_message
      use redistribute, only: scatter, gather
      use stella_time, only: code_dt
      use z_grid, only: nzgrid, ntubes
      use species, only: spec
      use parameters_physics, only: fphi
      use parameters_physics, only: radial_variation, full_flux_surface
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: rho_d_clamped
      use calculations_kxky, only: multiply_by_rho
      use velocity_grids, only: nvpa, nmu
      use velocity_grids, only: set_vpa_weights
      use geometry, only: bmag, dBdrho
      use stella_layouts, only: vmu_lo, kxkyz_lo
      use stella_layouts, only: is_idx, iky_idx, ikx_idx, iz_idx
      use dist_redistribute, only: kxkyz2vmu
      use arrays_dist_fn, only: gvmu, kperp2, dkperp2dr
      use arrays_fields, only: phi_corr_QN
      use g_tofrom_h, only: g_to_h
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phi, bpar
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gke_rhs
      real, dimension(:, :), intent(in out) :: time_collisions

      integer :: is, ikxkyz, imu, iv, ivmu, ikx, iky, iz, ia, it
      logical :: conservative_wgts
      real :: tfac, kfac

      complex, dimension(:), allocatable :: mucoll
      complex, dimension(:, :, :), allocatable :: coll
      complex, dimension(:, :, :, :, :), allocatable :: tmp_vmulo

      complex, dimension(:, :), allocatable :: g0k, g0x

      ia = 1

      if (full_flux_surface) then
         call mp_abort("collisions not currently supported for full_flux_surface=T.  Aborting.")
      end if

      if (proc0) call time_message(.false., time_collisions(:, 1), ' collisions')

      kfac = 0.0
      if (mu_operator) kfac = kfac + 0.5
      if (vpa_operator) kfac = kfac + 0.5

      allocate (tmp_vmulo(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      ! want exact conservation properties for collision operator
      conservative_wgts = .true.
      call set_vpa_weights(conservative_wgts)

      if (radial_variation) then
         allocate (g0k(naky, nakx))
         allocate (g0x(naky, nakx))
         !TODO (DSO) - could perhaps operator split the profile variation pieces off the main pieces, and so
         !             this portion of the code could just treat the terms that vary in x

         ! switch from g = <f> to h = f + Z*e*phi/T * F0
         tmp_vmulo = g
         call g_to_h(tmp_vmulo, phi, bpar, fphi, phi_corr_QN)

         !handle gyroviscous term
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  g0k = 0.5 * kfac * tmp_vmulo(:, :, iz, it, ivmu) * kperp2(:, :, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2
                  gke_rhs(:, :, iz, it, ivmu) = gke_rhs(:, :, iz, it, ivmu) - code_dt * spec(is)%vnew(is) * g0k

                  g0k = g0k * (dkperp2dr(:, :, ia, iz) - 2.0 * dBdrho(iz) / bmag(ia, iz) - spec(is)%tprim)
                  call multiply_by_rho(g0k)
                  gke_rhs(:, :, iz, it, ivmu) = gke_rhs(:, :, iz, it, ivmu) - code_dt * spec(is)%vnew(is) * g0k
               end do
            end do
         end do

         !handle the conservation terms
         if (momentum_conservation) call conserve_momentum_vmulo(tmp_vmulo, gke_rhs)
         if (energy_conservation) call conserve_energy_vmulo(tmp_vmulo, gke_rhs)

         !since Bessel functions do not appear under the velocity derivatives, these terms are one-point in x space
         ! and we can simply inverse Fourier transform
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call transform_kx2x_unpadded(tmp_vmulo(:, :, iz, it, ivmu), g0x)
                  tmp_vmulo(:, :, iz, it, ivmu) = g0x
               end do
            end do
         end do

         ! remap so that (vpa,mu) local
         if (proc0) call time_message(.false., time_collisions(:, 2), ' coll_redist')
         call scatter(kxkyz2vmu, tmp_vmulo, gvmu)
         if (proc0) call time_message(.false., time_collisions(:, 2), ' coll_redist')

         ! take vpa derivatives
         allocate (coll(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         allocate (mucoll(nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            if (vpa_operator) then
               !fix the temperature term
               tfac = (spec(is)%temp / spec(is)%temp_psi0) * (1.0 - rho_d_clamped(ikx) * spec(is)%tprim)
               do imu = 1, nmu
                  call vpa_differential_operator(tfac, gvmu(:, imu, ikxkyz), coll(:, imu, ikxkyz))
               end do
            else
               coll(:, :, ikxkyz) = 0.0
            end if

            if (mu_operator) then
               !fix the temperature/bmag term
               tfac = (spec(is)%temp / spec(is)%temp_psi0) &
                      * (1.0 - rho_d_clamped(ikx) * (spec(is)%tprim + dBdrho(iz) / bmag(ia, iz)))
               do iv = 1, nvpa
                  call mu_differential_operator(tfac, iz, ia, gvmu(iv, :, ikxkyz), mucoll)
                  coll(iv, :, ikxkyz) = coll(iv, :, ikxkyz) + mucoll
               end do
            end if
            gvmu(:, :, ikxkyz) = coll(:, :, ikxkyz)
         end do
         deallocate (coll, mucoll)

         ! remap so that (ky,kx,z,tube) local
         call gather(kxkyz2vmu, gvmu, tmp_vmulo)

         !don't forget to Fourier transform
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call transform_x2kx_unpadded(tmp_vmulo(:, :, iz, it, ivmu), g0k)
                  tmp_vmulo(:, :, iz, it, ivmu) = g0k
               end do
            end do
         end do

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            gke_rhs(:, :, :, :, ivmu) = gke_rhs(:, :, :, :, ivmu) + code_dt * spec(is)%vnew(is) * tmp_vmulo(:, :, :, :, ivmu)
         end do

         deallocate (g0k, g0x)
      else
         ! switch from g = <f> to h = f + Z*e*phi/T * F0
         tmp_vmulo = g
         call g_to_h(tmp_vmulo, phi, bpar, fphi)

         ! remap so that (vpa,mu) local
         if (proc0) call time_message(.false., time_collisions(:, 2), ' coll_redist')
         call scatter(kxkyz2vmu, tmp_vmulo, gvmu)
         if (proc0) call time_message(.false., time_collisions(:, 2), ' coll_redist')

         ia = 1

         ! take vpa derivatives
         allocate (coll(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
         allocate (mucoll(nmu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            if (vpa_operator) then
               do imu = 1, nmu
                  call vpa_differential_operator(1.0, gvmu(:, imu, ikxkyz), coll(:, imu, ikxkyz))
               end do
            else
               coll(:, :, ikxkyz) = 0.0
            end if
            if (mu_operator) then
               do iv = 1, nvpa
                  call mu_differential_operator(1.0, iz, ia, gvmu(iv, :, ikxkyz), mucoll)
                  coll(iv, :, ikxkyz) = coll(iv, :, ikxkyz) + mucoll
               end do
            end if
            if (momentum_conservation) call conserve_momentum(iky, ikx, iz, is, ikxkyz, gvmu(:, :, ikxkyz), coll(:, :, ikxkyz))
            if (energy_conservation) call conserve_energy(iz, is, ikxkyz, gvmu(:, :, ikxkyz), coll(:, :, ikxkyz))
            ! save memory by using gvmu and deallocating coll below
            ! before re-allocating tmp_vmulo
            gvmu(:, :, ikxkyz) = coll(:, :, ikxkyz) - 0.5 * kfac * kperp2(iky, ikx, ia, iz) * (spec(is)%smz / bmag(ia, iz))**2 * gvmu(:, :, ikxkyz)
         end do
         deallocate (coll, mucoll)

         ! remap so that (ky,kx,z,tube) local
         call gather(kxkyz2vmu, gvmu, tmp_vmulo)

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            gke_rhs(:, :, :, :, ivmu) = gke_rhs(:, :, :, :, ivmu) + code_dt * spec(is)%vnew(is) * tmp_vmulo(:, :, :, :, ivmu)
         end do
      end if

      deallocate (tmp_vmulo)

      ! reset to default integration wgts
      conservative_wgts = .false.
      call set_vpa_weights(conservative_wgts)

      if (proc0) call time_message(.false., time_collisions(:, 1), ' collisions')

   end subroutine advance_collisions_dougherty_explicit

   subroutine vpa_differential_operator(tfac, h, Dh)

      use velocity_grids, only: nvpa, vpa, dvpa

      implicit none

      real, intent(in) :: tfac
      complex, dimension(:), intent(in) :: h
      complex, dimension(:), intent(out) :: Dh

      integer :: iv

      ! use h = 0 at ghost cells beyond +/- vpa_max
      iv = 1
      Dh(iv) = (0.5 * h(iv + 1) * (tfac / dvpa + vpa(iv + 1)) - tfac * h(iv) / dvpa) / dvpa
      iv = nvpa
      Dh(iv) = (-tfac * h(iv) / dvpa + 0.5 * h(iv - 1) * (tfac / dvpa - vpa(iv - 1))) / dvpa
      do iv = 2, nvpa - 1
         Dh(iv) = (0.5 * h(iv + 1) * (tfac / dvpa + vpa(iv + 1)) - tfac * h(iv) / dvpa &
                   + 0.5 * h(iv - 1) * (tfac / dvpa - vpa(iv - 1))) / dvpa
      end do

   end subroutine vpa_differential_operator

   subroutine mu_differential_operator(tfac, iz, ia, h, Dh)

      use velocity_grids, only: nmu, dmu
      use velocity_grids, only: mu_cell, dmu_cell, wgts_mu_bare
      use velocity_grids, only: equally_spaced_mu_grid
      use geometry, only: bmag

      implicit none

      real, intent(in) :: tfac
      integer, intent(in) :: iz, ia
      complex, dimension(:), intent(in) :: h
      complex, dimension(:), intent(out) :: Dh

      integer :: imu
      real :: mm, m0, mp

      ! the following finite difference method is explained in init_mudiff_matrix

      imu = 1
      m0 = dmu_cell(imu) * (1.0 - tfac / (dmu(imu) * bmag(ia, iz))) / wgts_mu_bare(imu)
      mp = mu_cell(imu) * (tfac / (bmag(ia, iz) * dmu(imu)) + 1.0) / wgts_mu_bare(imu)
      Dh(imu) = m0 * h(imu) + mp * h(imu + 1)

      imu = nmu
      mm = mu_cell(imu - 1) * (tfac / (bmag(ia, iz) * dmu(imu - 1)) - 1.0) / wgts_mu_bare(imu)
      m0 = -mu_cell(imu - 1) * (tfac / (dmu(imu - 1) * bmag(ia, iz)) + 1.0) / wgts_mu_bare(imu)
      Dh(imu) = mm * h(imu - 1) + m0 * h(imu)

      do imu = 2, nmu - 1
         mm = mu_cell(imu - 1) * (tfac / (bmag(ia, iz) * dmu(imu - 1)) - 1.0) / wgts_mu_bare(imu)
         m0 = -(mu_cell(imu) / dmu(imu) + mu_cell(imu - 1) / dmu(imu - 1)) &
              * tfac / (wgts_mu_bare(imu) * bmag(ia, iz)) + dmu_cell(imu) / wgts_mu_bare(imu)
         mp = mu_cell(imu) * (tfac / (bmag(ia, iz) * dmu(imu)) + 1.0) / wgts_mu_bare(imu)
         Dh(imu) = mm * h(imu - 1) + m0 * h(imu) + mp * h(imu + 1)
      end do

   end subroutine mu_differential_operator

   subroutine conserve_momentum(iky, ikx, iz, is, ikxkyz, h, Ch)

      use species, only: spec
      use geometry, only: bmag
      use velocity_grids, only: integrate_vmu
      use velocity_grids, only: vpa, nvpa, nmu, vperp2
      use velocity_grids, only: maxwell_vpa, maxwell_mu
!     use velocity_grids, only: int_vpa2
      use arrays_dist_fn, only: kperp2
      use gyro_averages, only: aj0v, aj1v

      implicit none

      integer, intent(in) :: iky, ikx, iz, is, ikxkyz
      complex, dimension(:, :), intent(in) :: h
      complex, dimension(:, :), intent(in out) :: Ch

      complex, dimension(:, :), allocatable :: u_fac
      complex :: integral
      integer :: ia

      real :: norm

      allocate (u_fac(nvpa, nmu))

      ia = 1

!   norm = 1.0/int_vpa2(ia,iz,is)
      norm = 2.0

      if (vpa_operator) then
         u_fac = spread(aj0v(:, ikxkyz), 1, nvpa) * spread(vpa, 2, nmu)
         call integrate_vmu(u_fac * h, iz, integral)

         Ch = Ch + norm * u_fac * integral * spread(maxwell_mu(1, iz, :, is), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu)
      end if

      if (mu_operator) then
         u_fac = spread(vperp2(ia, iz, :) * aj1v(:, ikxkyz), 1, nvpa) * sqrt(kperp2(iky, ikx, ia, iz)) * spec(is)%smz / bmag(ia, iz)
         call integrate_vmu(u_fac * h, iz, integral)

         Ch = Ch + norm * u_fac * integral * spread(maxwell_mu(1, iz, :, is), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu)
      end if

      deallocate (u_fac)

   end subroutine conserve_momentum

   subroutine conserve_energy(iz, is, ikxkyz, h, Ch)

      use velocity_grids, only: integrate_vmu
      use velocity_grids, only: vpa, nvpa, nmu, vperp2
      use velocity_grids, only: maxwell_vpa, maxwell_mu
!   use velocity_grids, only: int_unit, int_vpa2, int_vperp2, int_vfrth
      use gyro_averages, only: aj0v

      implicit none

      integer, intent(in) :: iz, is, ikxkyz
      complex, dimension(:, :), intent(in) :: h
      complex, dimension(:, :), intent(in out) :: Ch

      complex, dimension(:, :), allocatable :: T_fac
      complex :: integral
      real :: norm

      integer :: ia

      allocate (T_fac(nvpa, nmu))

      ia = 1
      T_fac = 0.0
      if (vpa_operator) T_fac = spread(aj0v(:, ikxkyz), 1, nvpa) * (spread(vpa**2, 2, nmu)) - 0.5
      if (mu_operator) T_fac = T_fac + spread(vperp2(ia, iz, :), 1, nvpa) - 1.0
      call integrate_vmu(T_fac * h, iz, integral)

!     norm = 1.0/(int_vfrth(ia,iz,is) - (int_vperp2(ia,iz,is) + int_vpa2(ia,iz,is))**2/int_unit(ia,iz,is))
      norm = 2.0 / 3.0

      Ch = Ch + 2.0 * norm * T_fac * integral * spread(maxwell_mu(ia, iz, :, is), 1, nvpa) * spread(maxwell_vpa(:, is), 2, nmu)

      deallocate (T_fac)

   end subroutine conserve_energy

   subroutine conserve_momentum_vmulo(h, gke_rhs)

      use mp, only: sum_allreduce
      use stella_time, only: code_dt
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, iv_idx, is_idx
      use species, only: spec
      use parameters_physics, only: radial_variation
      use geometry, only: bmag, dBdrho
      use parameters_kxky_grid, only: nakx, naky
      use calculations_kxky, only: multiply_by_rho
      use z_grid, only: nzgrid, ntubes
      use velocity_grids, only: integrate_species, mu, vpa, vperp2
      use velocity_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_dist_fn, only: kperp2, dkperp2dr
      use gyro_averages, only: gyro_average, gyro_average_j1, aj0x, aj1x

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: h
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gke_rhs

      complex, dimension(:, :), allocatable :: g0k, g1k
      complex, dimension(:, :, :), allocatable :: gyro_g
      complex, dimension(:, :, :, :), allocatable :: field1, field2
      real :: prefac, energy
      integer :: it, iz, ivmu, imu, iv, ia, is

      ia = 1

      allocate (g0k(naky, nakx))
      allocate (g1k(naky, nakx))
      allocate (gyro_g(naky, nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      allocate (field1(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (field2(naky, nakx, -nzgrid:nzgrid, ntubes))

      !component from vpa
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)

               call gyro_average(h(:, :, iz, it, ivmu), iz, ivmu, gyro_g(:, :, ivmu))
               gyro_g(:, :, ivmu) = gyro_g(:, :, ivmu) * vpa(iv)
               g0k = 0.0
               if (radial_variation) then
                  g0k(:, :) = gyro_g(:, :, ivmu) &
                              * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                 * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                 * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                 + dBdrho(iz) / bmag(ia, iz))

                  call multiply_by_rho(g0k)
               end if
               gyro_g(:, :, ivmu) = gyro_g(:, :, ivmu) + g0k

            end do
            call integrate_species(gyro_g, iz, spec%dens_psi0 * spec%temp_psi0, field1(:, :, iz, it), reduce_in=.false.)
         end do
      end do
      call sum_allreduce(field1)

      !component from vperp
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)

               !component from vpa
               call gyro_average_j1(h(:, :, iz, it, ivmu), iz, ivmu, gyro_g(:, :, ivmu))
               gyro_g(:, :, ivmu) = gyro_g(:, :, ivmu) * vperp2(ia, iz, imu) * sqrt(kperp2(:, :, ia, iz)) * spec(is)%smz_psi0 / bmag(ia, iz)
               g0k = 0.0
               if (radial_variation) then
                  g0k = gyro_g(:, :, ivmu) * (dBdrho(iz) / bmag(ia, iz) + 0.5 * dkperp2dr(:, :, ia, iz)) &
                        + h(:, :, iz, it, ivmu) * (0.5 * aj0x(:, :, iz, ivmu) - aj1x(:, :, iz, ivmu)) &
                        * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) * vperp2(ia, iz, imu) &
                        * sqrt(kperp2(:, :, ia, iz)) * spec(is)%smz_psi0 / bmag(ia, iz)

                  call multiply_by_rho(g0k)
               end if
               gyro_g(:, :, ivmu) = gyro_g(:, :, ivmu) + g0k

            end do
            call integrate_species(gyro_g, iz, spec%dens_psi0 * spec%temp_psi0, field2(:, :, iz, it), reduce_in=.false.)

         end do
      end do
      call sum_allreduce(field2)
      deallocate (gyro_g)

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)

               prefac = 2.0 / (spec(is)%dens * spec(is)%temp) * code_dt * spec(is)%vnew(is) &
                        * maxwell_mu(ia, iz, imu, is) * maxwell_vpa(iv, is) * maxwell_fac(is)

               g0k = aj0x(:, :, iz, ivmu) * vpa(iv) * field1(:, :, iz, it) &
                     + aj1x(:, :, iz, ivmu) * vperp2(ia, iz, imu) * field2(:, :, iz, it) &
                     * spec(is)%smz_psi0 * sqrt(kperp2(:, :, ia, iz)) / bmag(ia, iz)

               gke_rhs(:, :, iz, it, ivmu) = gke_rhs(:, :, iz, it, ivmu) + prefac * g0k

               if (radial_variation) then
                  energy = (vpa(iv)**2 + vperp2(ia, iz, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)

                  g1k = field1(:, :, iz, it) * vpa(iv) * (-0.5 * aj1x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                                          * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                                          * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))) &
                        + field2(:, :, iz, it) * spec(is)%smz_psi0 * vperp2(ia, iz, imu) / bmag(ia, iz) * sqrt(kperp2(:, :, ia, iz)) &
                        * (0.5 * aj1x(:, :, iz, ivmu) * dkperp2dr(:, :, ia, iz) + (0.5 * aj0x(:, :, iz, ivmu) - aj1x(:, :, iz, ivmu)) &
                           * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)))
                  g1k = g1k + g0k * (spec(is)%tprim * (energy - 2.5) + 2 * mu(imu) * dBdrho(iz))

                  call multiply_by_rho(g1k)

                  gke_rhs(:, :, iz, it, ivmu) = gke_rhs(:, :, iz, it, ivmu) + prefac * g1k
               end if
            end do
         end do
      end do

      deallocate (g0k, g1k)
      deallocate (field1, field2)

   end subroutine conserve_momentum_vmulo

   subroutine conserve_energy_vmulo(h, gke_rhs)

      use mp, only: sum_allreduce
      use stella_time, only: code_dt
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: imu_idx, iv_idx, is_idx
      use species, only: spec
      use parameters_physics, only: radial_variation
      use geometry, only: bmag, dBdrho
      use parameters_kxky_grid, only: nakx, naky
      use calculations_kxky, only: multiply_by_rho
      use z_grid, only: nzgrid, ntubes
      use velocity_grids, only: integrate_species
      use velocity_grids, only: mu, vpa, nmu, vperp2
      use velocity_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use arrays_dist_fn, only: kperp2, dkperp2dr
      use gyro_averages, only: gyro_average, gyro_average_j1, aj0x, aj1x

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: h
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gke_rhs

      complex, dimension(:, :), allocatable :: g0k, g1k
      complex, dimension(:, :, :), allocatable :: gyro_g
      complex, dimension(:, :, :, :), allocatable :: field
      real :: prefac, energy
      integer :: it, iz, ivmu, imu, iv, ia, is

      ia = 1

      allocate (g0k(naky, nakx))
      allocate (g1k(naky, nakx))
      allocate (gyro_g(naky, nakx, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      allocate (field(naky, nakx, -nzgrid:nzgrid, ntubes))

      !component from vpa
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)

               call gyro_average(h(:, :, iz, it, ivmu), iz, ivmu, gyro_g(:, :, ivmu))
               gyro_g(:, :, ivmu) = gyro_g(:, :, ivmu) * (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5 * spec(is)%temp / spec(is)%temp_psi0)
               g0k = 0.0
               if (radial_variation) then
                  g0k(:, :) = gyro_g(:, :, ivmu) &
                              * (-0.5 * aj1x(:, :, iz, ivmu) / aj0x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                 * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                 * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz)) &
                                 + dBdrho(iz) / bmag(ia, iz)) &
                              + h(:, :, iz, it, ivmu) * aj0x(:, :, iz, ivmu) * (vperp2(ia, iz, imu) * dBdrho(iz) / bmag(ia, iz) &
                                                                                + 1.5 * spec(is)%tprim)
                  call multiply_by_rho(g0k)
               end if
               gyro_g(:, :, ivmu) = gyro_g(:, :, ivmu) + g0k

            end do
            call integrate_species(gyro_g, iz, spec%dens_psi0 * spec%temp_psi0**2, field(:, :, iz, it), reduce_in=.false.)
         end do
      end do
      call sum_allreduce(field)

      deallocate (gyro_g)

      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               is = is_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               iv = iv_idx(vmu_lo, ivmu)

               prefac = 4.0 / (3.0 * spec(is)%dens * spec(is)%temp**2) * code_dt * spec(is)%vnew(is) &
                        * maxwell_mu(ia, iz, imu, is) * maxwell_vpa(iv, is) * maxwell_fac(is)

               g0k = aj0x(:, :, iz, ivmu) * field(:, :, iz, it) &
                     * (vpa(iv)**2 + vperp2(ia, iz, imu) - 1.5 * spec(is)%temp / spec(is)%temp_psi0)

               gke_rhs(:, :, iz, it, ivmu) = gke_rhs(:, :, iz, it, ivmu) + prefac * g0k

               if (radial_variation) then
                  energy = (vpa(iv)**2 + vperp2(ia, iz, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)

                  g1k = field(:, :, iz, it) * (energy - 1.5) * (-0.5 * aj1x(:, :, iz, ivmu) * (spec(is)%smz)**2 &
                                                                * (kperp2(:, :, ia, iz) * vperp2(ia, iz, imu) / bmag(ia, iz)**2) &
                                                                * (dkperp2dr(:, :, ia, iz) - dBdrho(iz) / bmag(ia, iz))) &
                        + field(:, :, iz, it) * aj0x(:, :, iz, ivmu) * (vperp2(ia, iz, imu) * dBdrho(iz) / bmag(ia, iz) + 1.5 * spec(is)%tprim)

                  g1k = g1k + g0k * (spec(is)%tprim * (energy - 3.5) + 2 * mu(imu) * dBdrho(iz))

                  call multiply_by_rho(g1k)

                  gke_rhs(:, :, iz, it, ivmu) = gke_rhs(:, :, iz, it, ivmu) + prefac * g1k
               end if
            end do
         end do
      end do

      deallocate (g0k, g1k)
      deallocate (field)

   end subroutine conserve_energy_vmulo

   subroutine advance_collisions_dougherty_implicit(phi, apar, bpar)

      use z_grid, only: nzgrid
      use arrays_dist_fn, only: gvmu

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar

      if (vpa_operator) call advance_vpadiff_implicit(phi, apar, bpar, gvmu)
      if (mu_operator) call advance_mudiff_implicit(phi, apar, bpar, gvmu)

   end subroutine advance_collisions_dougherty_implicit

   subroutine advance_vpadiff_implicit(phi, apar, bpar, g)

      use mp, only: sum_allreduce
      use finite_differences, only: tridag
      use linear_solve, only: lu_back_substitution
      use stella_time, only: code_dt
      use parameters_physics, only: fphi
      use species, only: nspec, spec, has_electron_species
      use z_grid, only: nzgrid, ntubes
      use velocity_grids, only: nmu, nvpa
      use velocity_grids, only: maxwell_vpa, maxwell_mu, vpa, vperp2
      use velocity_grids, only: set_vpa_weights
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: zonal_mode
      use geometry, only: dl_over_b
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use g_tofrom_h, only: g_to_h
      use gyro_averages, only: aj0v
      use fields_fluxtube, only: get_fields_fluxtube
      use arrays_fields, only: efac, gamtot_h
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g

      integer :: ikxkyz, iky, ikx, iz, it, is, ia
      integer :: imu
      integer :: idx
      real, dimension(:, :), allocatable :: tmp
      complex, dimension(:, :, :, :, :), allocatable :: flds
      complex, dimension(:), allocatable :: tmp2
      complex, dimension(:, :, :), allocatable :: flds_zf, g_in

      ia = 1

      ! store input g for use later, as g will be overwritten below
      allocate (g_in(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      g_in = g

      ! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
      ! g = g^{***}.  tridag below inverts above equation to get h_inh^{n+1}
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            call tridag(1, aa_vpa(:, is), bb_vpa(:, ikxkyz), cc_vpa(:, is), g(:, imu, ikxkyz))
         end do
      end do

      allocate (flds(naky, nakx, -nzgrid:nzgrid, ntubes, nresponse_vpa))
      allocate (flds_zf(nakx, ntubes, nresponse_vpa)); flds_zf = 0.

      ! need to obtain phi^{n+1} and conservation terms using response matrix approach
      ! first get phi_inh^{n+1}
      call get_fields_fluxtube(g, phi, apar, bpar, dist='h', skip_fsa=.true.)
      flds(:, :, :, :, 1) = phi

      idx = 2
      ! get upar_inh^{n+1}
      if (momentum_conservation) then
         call get_upar(g, flds(:, :, :, :, idx:idx + nspec - 1))
         idx = idx + nspec
      end if

      ! get temp_inh^{n+1}
      if (energy_conservation) call get_temp(g, flds(:, :, :, :, idx:idx + nspec - 1))

      phi = 0.0
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         ! all is indices inside ikxkyz super-index have same info
         ! no need to compute multiple times
         is = is_idx(kxkyz_lo, ikxkyz); if (is /= 1) cycle
         call lu_back_substitution(vpadiff_response(:, :, ikxkyz), vpadiff_idx(:, ikxkyz), &
                                   flds(iky, ikx, iz, it, :))
         if (.not. has_electron_species(spec) .and. zonal_mode(iky) &
             .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            flds_zf(ikx, it, :) = flds_zf(ikx, it, :) + dl_over_b(ia, iz) * flds(iky, ikx, iz, it, :)
         end if
         phi(iky, ikx, iz, it) = flds(iky, ikx, iz, it, 1)
      end do

      if (.not. has_electron_species(spec) .and. zonal_mode(1) &
          .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         !complete flux-surface average
         call sum_allreduce(flds_zf)
         do it = 1, ntubes
            do ikx = 1, nakx
               call lu_back_substitution(vpadiff_zf_response(:, :, ikx), vpadiff_zf_idx(:, ikx), &
                                         flds_zf(ikx, it, :))
               !multiply by Q, which has a single non-zero component
               flds_zf(ikx, it, 1) = (efac / gamtot_h) * flds_zf(ikx, it, 1)
               flds_zf(ikx, it, 2:) = 0.
            end do
         end do

         allocate (tmp2(nresponse_vpa))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            if (iky /= 1 .or. is /= 1) cycle

            tmp2 = flds_zf(ikx, it, :)
            call lu_back_substitution(vpadiff_response(:, :, ikxkyz), vpadiff_idx(:, ikxkyz), &
                                      tmp2)

            phi(1, ikx, iz, it) = phi(1, ikx, iz, it) + tmp2(1)
         end do
         deallocate (tmp2)

      end if

      call sum_allreduce(phi)

      g = g_in

      ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + 2*dt*nu*J0*F0*(vpa*upar+(v^2-3/2)*temp)
      ! first two terms added via g_to_h subroutine
      call g_to_h(g, phi, bpar, fphi)

      allocate (tmp(nvpa, nmu))

      if (momentum_conservation .or. energy_conservation) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            tmp = 2.0 * code_dt * spec(is)%vnew(is) &
                  * spread(maxwell_vpa(:, is), 2, nmu) * spread(aj0v(:, ikxkyz) * maxwell_mu(1, iz, :, is), 1, nvpa)
            if (momentum_conservation) &
               g(:, :, ikxkyz) = g(:, :, ikxkyz) + tmp * spread(vpa * flds(iky, ikx, iz, it, is + 1), 2, nmu)
            if (energy_conservation) &
               g(:, :, ikxkyz) = g(:, :, ikxkyz) &
                                 + tmp * (spread(vpa**2, 2, nmu) + spread(vperp2(1, iz, :), 1, nvpa) - 1.5) * flds(iky, ikx, iz, it, idx + is - 1)
         end do
      end if

      deallocate (tmp, flds, flds_zf)

      ! now invert system to get h^{n+1}
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         is = is_idx(kxkyz_lo, ikxkyz)
         do imu = 1, nmu
            call tridag(1, aa_vpa(:, is), bb_vpa(:, ikxkyz), cc_vpa(:, is), g(:, imu, ikxkyz))
         end do
      end do

      ! now get g^{n+1} from h^{n+1} and phi^{n+1}
      call g_to_h(g, phi, bpar, -fphi)

      deallocate (g_in)

   end subroutine advance_vpadiff_implicit

   subroutine advance_mudiff_implicit(phi, apar, bpar, g)

      use mp, only: sum_allreduce
      use finite_differences, only: tridag
      use linear_solve, only: lu_back_substitution
      use stella_time, only: code_dt
      use parameters_physics, only: fphi
      use species, only: nspec, spec, has_electron_species
      use z_grid, only: nzgrid, ntubes
      use velocity_grids, only: nmu, nvpa
      use velocity_grids, only: maxwell_vpa, maxwell_mu, vpa, vperp2
      use velocity_grids, only: set_vpa_weights
      use parameters_kxky_grid, only: naky, nakx
      use grids_kxky, only: zonal_mode
      use stella_layouts, only: kxkyz_lo
      use stella_layouts, only: iky_idx, ikx_idx, iz_idx, it_idx, is_idx
      use arrays_dist_fn, only: kperp2
      use gyro_averages, only: aj0v, aj1v
      use g_tofrom_h, only: g_to_h
      use fields_fluxtube, only: get_fields_fluxtube
      use arrays_fields, only: efac, gamtot_h
      use geometry, only: bmag, dl_over_b
      use parameters_physics, only: adiabatic_option_switch
      use parameters_physics, only: adiabatic_option_fieldlineavg

      ! TMP FOR TESTING
!    use velocity_grids, only: mu

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar
      complex, dimension(:, :, kxkyz_lo%llim_proc:), intent(in out) :: g

      integer :: ikxkyz, iky, ikx, iz, it, is, ia
      integer :: iv
      integer :: idx

      ! TMP FOR TESTING
!    integer :: imu

      real, dimension(:, :), allocatable :: tmp
      complex, dimension(:), allocatable :: tmp2
      complex, dimension(:, :, :, :, :), allocatable :: flds
      complex, dimension(:, :, :), allocatable :: flds_zf, g_in

      ia = 1

      ! store input g for use later, as g will be overwritten below
      allocate (g_in(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      g_in = g

      ! since backwards difference in time, (I-dt*D)h_inh^{n+1} = g^{***}
      ! g = g^{***}.  tridag below inverts above equation to get h_inh^{n+1}
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         ! TMP FOR TESTING
!       do imu = 1, nmu
!          g(:,imu,ikxkyz) = maxwell_vpa*maxwell_mu(1,iz,imu)
!       end do
         do iv = 1, nvpa
            call tridag(1, aa_mu(iz, :, is), bb_mu(:, ikxkyz), cc_mu(iz, :, is), g(iv, :, ikxkyz))
         end do
         ! TMP FOR TESTING
!       iv = nvpa/2
!       do imu = 1, nmu
!          write (*,*) 'ggg', mu(imu), real(g(iv,imu,ikxkyz)), aimag(g(iv,imu,ikxkyz)), maxwell_vpa(iv,is)*maxwell_mu(1,iz,imu)
!       end do
      end do

      allocate (flds(naky, nakx, -nzgrid:nzgrid, ntubes, nresponse_mu))
      allocate (flds_zf(nakx, ntubes, nresponse_mu)); flds_zf = 0.

      ! need to obtain phi^{n+1} and conservation terms using response matrix approach
      ! first get phi_inh^{n+1}
      call get_fields_fluxtube(g, phi, apar, bpar, dist='h', skip_fsa=.true.)
      flds(:, :, :, :, 1) = phi

      idx = 2
      ! get upar_inh^{n+1}
      if (momentum_conservation) then
         call get_uperp(g, flds(:, :, :, :, idx:idx + nspec - 1))
         idx = idx + nspec
      end if

      ! get temp_inh^{n+1}
      if (energy_conservation) call get_temp_mu(g, flds(:, :, :, :, idx:idx + nspec - 1))

      phi = 0.0
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iky = iky_idx(kxkyz_lo, ikxkyz)
         ikx = ikx_idx(kxkyz_lo, ikxkyz)
         iz = iz_idx(kxkyz_lo, ikxkyz)
         it = it_idx(kxkyz_lo, ikxkyz)
         ! all is indices inside ikxkyz super-index have same info
         ! no need to compute multiple times
         is = is_idx(kxkyz_lo, ikxkyz); if (is /= 1) cycle
         call lu_back_substitution(mudiff_response(:, :, ikxkyz), mudiff_idx(:, ikxkyz), &
                                   flds(iky, ikx, iz, it, :))
         if (.not. has_electron_species(spec) .and. zonal_mode(iky) &
             .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            flds_zf(ikx, it, :) = flds_zf(ikx, it, :) + dl_over_b(ia, iz) * flds(iky, ikx, iz, it, :)
         end if
         phi(iky, ikx, iz, it) = flds(iky, ikx, iz, it, 1)
      end do

      if (.not. has_electron_species(spec) .and. zonal_mode(1) &
          .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         !complete flux-surface average
         call sum_allreduce(flds_zf)
         do it = 1, ntubes
            do ikx = 1, nakx
               call lu_back_substitution(mudiff_zf_response(:, :, ikx), mudiff_zf_idx(:, ikx), &
                                         flds_zf(ikx, it, :))
               !multiply by Q, which has a single non-zero component
               flds_zf(ikx, it, 1) = (efac / gamtot_h) * flds_zf(ikx, it, 1)
               flds_zf(ikx, it, 2:) = 0.
            end do
         end do

         allocate (tmp2(nresponse_mu))
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)

            if (iky /= 1 .or. is /= 1) cycle

            tmp2 = flds_zf(ikx, it, :)
            call lu_back_substitution(mudiff_response(:, :, ikxkyz), mudiff_idx(:, ikxkyz), &
                                      tmp2)

            phi(1, ikx, iz, it) = phi(1, ikx, iz, it) + tmp2(1)
         end do
         deallocate (tmp2)

      end if

      call sum_allreduce(phi)

      g = g_in

      ! RHS is g^{***} + Ze/T*<phi^{n+1}>*F0 + ...
      ! first two terms added via g_to_h subroutine
      call g_to_h(g, phi, bpar, fphi)

      allocate (tmp(nvpa, nmu))

      if (momentum_conservation .or. energy_conservation) then
         do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
            iky = iky_idx(kxkyz_lo, ikxkyz)
            ikx = ikx_idx(kxkyz_lo, ikxkyz)
            iz = iz_idx(kxkyz_lo, ikxkyz)
            it = it_idx(kxkyz_lo, ikxkyz)
            is = is_idx(kxkyz_lo, ikxkyz)
            tmp = 2.0 * code_dt * spec(is)%vnew(is) &
                  * spread(maxwell_vpa(:, is), 2, nmu) * spread(maxwell_mu(1, iz, :, is), 1, nvpa)
            if (momentum_conservation) &
               g(:, :, ikxkyz) = g(:, :, ikxkyz) + tmp * kperp2(iky, ikx, ia, iz) &
                                 * spread(vperp2(ia, iz, :) * aj1v(:, ikxkyz), 1, nvpa) * (spec(is)%smz / bmag(ia, iz))**2 &
                                 * flds(iky, ikx, iz, it, is + 1)
            if (energy_conservation) &
               g(:, :, ikxkyz) = g(:, :, ikxkyz) &
                                 + flds(iky, ikx, iz, it, idx + is - 1) * tmp * spread(aj0v(:, ikxkyz), 1, nvpa) &
                                 * (spread(vpa**2, 2, nmu) + spread(vperp2(1, iz, :), 1, nvpa) - 1.5)
         end do
      end if

      deallocate (tmp, flds, flds_zf)

      ! now invert system to get h^{n+1}
      do ikxkyz = kxkyz_lo%llim_proc, kxkyz_lo%ulim_proc
         iz = iz_idx(kxkyz_lo, ikxkyz)
         is = is_idx(kxkyz_lo, ikxkyz)
         do iv = 1, nvpa
            call tridag(1, aa_mu(iz, :, is), bb_mu(:, ikxkyz), cc_mu(iz, :, is), g(iv, :, ikxkyz))
         end do
      end do

      ! now get g^{n+1} from h^{n+1} and phi^{n+1}
      call g_to_h(g, phi, bpar, -fphi)

      deallocate (g_in)

   end subroutine advance_mudiff_implicit

end module coll_dougherty
