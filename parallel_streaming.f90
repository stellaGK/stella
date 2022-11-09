module parallel_streaming

   implicit none

   public :: init_parallel_streaming, finish_parallel_streaming
   public :: advance_parallel_streaming_explicit
   public :: advance_parallel_streaming_implicit
   public :: add_parallel_streaming_radial_variation
   public :: stream_tridiagonal_solve
   public :: get_gke_rhs
   public :: sweep_zed_zonal
   public :: parallel_streaming_initialized
   public :: stream, stream_c, stream_sign
   public :: time_parallel_streaming
   public :: stream_rad_var1
   public :: stream_rad_var2

   private

   interface center_zed
      module procedure center_zed_segment_real
      module procedure center_zed_extended
   end interface

   logical :: parallel_streaming_initialized = .false.

   integer, dimension(:), allocatable :: stream_sign
   real, dimension(:, :, :, :), allocatable :: stream
   real, dimension(:, :, :), allocatable :: stream_c
   real, dimension(:, :, :), allocatable :: stream_rad_var1
   real, dimension(:, :, :), allocatable :: stream_rad_var2
   real, dimension(:, :), allocatable :: stream_tri_a1, stream_tri_a2
   real, dimension(:, :), allocatable :: stream_tri_b1, stream_tri_b2
   real, dimension(:, :), allocatable :: stream_tri_c1, stream_tri_c2
   real, dimension(:, :), allocatable :: drift_tri_a1
   real, dimension(:, :), allocatable :: drift_tri_b1
   real, dimension(:, :), allocatable :: drift_tri_c1
   real, dimension(:, :), allocatable :: gradpar_c

   real, dimension(2, 3) :: time_parallel_streaming = 0.

contains

   subroutine init_parallel_streaming

      use finite_differences, only: fd3pt
      use stella_time, only: code_dt
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use species, only: spec, nspec, pfac
      use vpamu_grids, only: nvpa, nvpa
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: vperp2, vpa, mu
      use kt_grids, only: nalpha
      use zgrid, only: nzgrid, nztot
      use stella_geometry, only: gradpar, dgradpardrho, dBdrho, gfac, b_dot_grad_z
      use run_parameters, only: stream_implicit, driftkinetic_implicit
      use physics_flags, only: include_parallel_streaming, radial_variation

      implicit none

      integer :: iv, imu, is, ivmu
      integer :: ia, iz

      real, dimension(:), allocatable :: energy

      if (parallel_streaming_initialized) return
      parallel_streaming_initialized = .true.

      if (.not. allocated(stream)) allocate (stream(nalpha, -nzgrid:nzgrid, nvpa, nspec)); stream = 0.
      if (.not. allocated(stream_sign)) allocate (stream_sign(nvpa)); stream_sign = 0

      ! sign of stream corresponds to appearing on RHS of GK equation
      ! i.e., this is the factor multiplying dg/dz on RHS of equation
      if (include_parallel_streaming) then
         do iv = 1, nvpa
            do iz = -nzgrid, nzgrid
               do ia = 1, nalpha
                  stream(ia, iz, iv, :) = -code_dt * b_dot_grad_z(ia, iz) * vpa(iv) * spec%stm_psi0
               end do
            end do
         end do
      else
         stream = 0.0
      end if

      if (radial_variation) then
         allocate (energy(-nzgrid:nzgrid))

         if (.not. allocated(stream_rad_var1)) then
            allocate (stream_rad_var1(-nzgrid:nzgrid, nvpa, nspec))
         end if
         if (.not. allocated(stream_rad_var2)) then
            allocate (stream_rad_var2(nalpha, -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
            stream_rad_var2 = 0.0
         end if
         ia = 1
         stream_rad_var1 = -code_dt * spread(spread(spec%stm_psi0, 1, nztot), 2, nvpa) &
                           * gfac * spread(spread(vpa, 1, nztot) * spread(dgradpardrho, 2, nvpa), 3, nspec)
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            energy = (vpa(iv)**2 + vperp2(ia, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
            stream_rad_var2(ia, :, ivmu) = &
               +code_dt * spec(is)%stm_psi0 * vpa(iv) * gradpar &
               * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_mu(ia, :, imu, is) * maxwell_fac(is) &
               * (pfac * (spec(is)%fprim + spec(is)%tprim * (energy - 2.5)) &
                  + gfac * 2 * mu(imu) * dBdrho)
         end do
         deallocate (energy)
      end if

      !> stream_sign set to +/- 1 depending on the sign of the parallel streaming term.
      !> NB: stream_sign = -1 corresponds to positive advection velocity
      !> only need to consider ia=1, iz=0 and is=1 because alpha, z and species dependences
      !> do not lead to change in sign of the streaming pre-factor
      do iv = 1, nvpa
         stream_sign(iv) = int(sign(1.0, stream(1, 0, iv, 1)))
      end do

      if (stream_implicit .or. driftkinetic_implicit) then
         call init_invert_stream_operator
         if (.not. allocated(stream_c)) allocate (stream_c(-nzgrid:nzgrid, nvpa, nspec))
         stream_c = stream(1, :, :, :)
         do is = 1, nspec
            do iv = 1, nvpa
               call center_zed(iv, stream_c(:, iv, is))
            end do
         end do
         if (.not. allocated(gradpar_c)) allocate (gradpar_c(-nzgrid:nzgrid, -1:1))
         gradpar_c = spread(gradpar, 2, 3)
         !> get gradpar centred in zed for negative vpa (affects upwinding)
         call center_zed(1, gradpar_c(:, -stream_sign(1)))
         !> get gradpar centred in zed for positive vpa (affects upwinding)
         call center_zed(nvpa, gradpar_c(:, -stream_sign(nvpa)))
         stream = spread(stream_c, 1, nalpha)
      end if

   end subroutine init_parallel_streaming

   subroutine init_invert_stream_operator

      use zgrid, only: delzed
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: nsegments
      use run_parameters, only: zed_upwind, time_upwind, stream_drifts_implicit

      implicit none

      integer :: nz, nseg_max

      nz = maxval(iz_up - iz_low)
      nseg_max = maxval(nsegments)

      if (.not. allocated(stream_tri_a1)) then
         allocate (stream_tri_a1(nz * nseg_max + 1, -1:1)); stream_tri_a1 = 0.
         allocate (stream_tri_a2(nz * nseg_max + 1, -1:1)); stream_tri_a2 = 0.
         allocate (stream_tri_b1(nz * nseg_max + 1, -1:1)); stream_tri_b1 = 1.
         allocate (stream_tri_b2(nz * nseg_max + 1, -1:1)); stream_tri_b2 = 0.
         allocate (stream_tri_c1(nz * nseg_max + 1, -1:1)); stream_tri_c1 = 0.
         allocate (stream_tri_c2(nz * nseg_max + 1, -1:1)); stream_tri_c2 = 0.
      end if

      ! corresponds to sign of stream term positive on RHS of equation
      ! i.e., negative parallel advection speed
      ! NB: assumes equal spacing in zed
      stream_tri_b1(:, 1) = 0.5 * (1.0 + zed_upwind)
      stream_tri_b2(:, 1) = -1.0 / delzed(0)
      stream_tri_c1(:nz * nseg_max, 1) = 0.5 * (1.0 - zed_upwind)
      stream_tri_c2(:nz * nseg_max, 1) = 1.0 / delzed(0)

      ! corresponds to sign of stream term negative on RHS of equation
      ! NB: assumes equal spacing in zed
      stream_tri_b1(:, -1) = 0.5 * (1.0 + zed_upwind)
      stream_tri_b2(:, -1) = 1.0 / delzed(0)
      stream_tri_a1(2:, -1) = 0.5 * (1.0 - zed_upwind)
      stream_tri_a2(2:, -1) = -1.0 / delzed(0)

      stream_tri_a2 = 0.5 * (1.0 + time_upwind) * stream_tri_a2
      stream_tri_b2 = 0.5 * (1.0 + time_upwind) * stream_tri_b2
      stream_tri_c2 = 0.5 * (1.0 + time_upwind) * stream_tri_c2

      if (stream_drifts_implicit) then
         if (.not. allocated(drift_tri_a1)) then
            allocate (drift_tri_a1(nz * nseg_max + 1, -1:1)); drift_tri_a1 = 0.
            allocate (drift_tri_b1(nz * nseg_max + 1, -1:1)); drift_tri_b1 = 0.
            allocate (drift_tri_c1(nz * nseg_max + 1, -1:1)); drift_tri_c1 = 0.
         end if

         ! corresponds to sign of stream term positive on RHS of equation
         ! i.e., negative parallel advection speed
         drift_tri_b1(:, 1) = 0.5 * (1.0 + zed_upwind)
         drift_tri_c1(:nz * nseg_max, 1) = 0.5 * (1.0 - zed_upwind)

         ! corresponds to sign of stream term negative on RHS of equation
         drift_tri_b1(:, -1) = 0.5 * (1.0 + zed_upwind)
         drift_tri_a1(2:, -1) = 0.5 * (1.0 - zed_upwind)

         drift_tri_a1 = 0.5 * (1.0 + time_upwind) * drift_tri_a1
         drift_tri_b1 = 0.5 * (1.0 + time_upwind) * drift_tri_b1
         drift_tri_c1 = 0.5 * (1.0 + time_upwind) * drift_tri_c1
      end if

   end subroutine init_invert_stream_operator

   subroutine advance_parallel_streaming_explicit(g, gout)

      use mp, only: proc0, mp_abort
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use job_manage, only: time_message
      use stella_transforms, only: transform_ky2y
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, naky_all, nakx, ikx_max, ny
      use kt_grids, only: swap_kxky
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use species, only: spec
      use physics_flags, only: full_flux_surface
      use gyro_averages, only: gyro_average
      use run_parameters, only: driftkinetic_implicit
      use run_parameters, only: fapar, fbpar
      use fields_arrays, only: phi, apar, bpar
      use fields, only: get_gyroaverage_chi

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      integer :: ivmu, iv, imu, is, ia, iz, it
      complex, dimension(:, :, :, :), allocatable :: g0, dgchi_dz
      complex, dimension(:, :, :, :), allocatable :: g0y, g1y
      complex, dimension(:, :), allocatable :: g0_swap

      !> if flux tube simulation parallel streaming stays in ky,kx,z space with ky,kx,z local
      !> if full flux surface (flux annulus), will need to calculate in y space

      !> start the timer for the parallel streaming part of the time advance
      if (proc0) call time_message(.false., time_parallel_streaming(:, 1), ' Stream advance')

      !> allocate arrays needed for intermmediate calculations
      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (dgchi_dz(naky, nakx, -nzgrid:nzgrid, ntubes))
      !> if simulating a full flux surface, will also need version of the above arrays
      !> that is Fourier transformed to y-space
      if (full_flux_surface) then
         allocate (g0_swap(naky_all, ikx_max))
         allocate (g0y(ny, ikx_max, -nzgrid:nzgrid, ntubes))
         allocate (g1y(ny, ikx_max, -nzgrid:nzgrid, ntubes))
      end if

      if (driftkinetic_implicit .and. ((fapar > epsilon(0.0)) .or. (fbpar > epsilon(0.0)))) then
         call mp_abort("driftkinetic_implicit not currently compatible with fapar>0 or fbpar>0. Aborting.")
      end if

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         !> get (iv,imu,is) indices corresponding to ivmu super-index
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! !!! Bob: commented out for h source term
         ! !> obtain <phi> (or <phi>-phi if driftkinetic_implicit=T)
         ! call get_gyroaverage_chi(ivmu, phi, apar, bpar, g0(:, :, :, :))
         ! if (driftkinetic_implicit) g0(:, :, :, :) = g0(:, :, :, :) - phi

         ! !> get d<chi>/dz, with z the parallel coordinate and store in dgchi_dz
         ! !> note that this should be a centered difference to avoid numerical
         ! !> unpleasantness to do with inexact cancellations in later velocity integration
         ! !> see appendix of the stella JCP 2019 for details
         ! call get_dgdz_centered(g0, ivmu, dgchi_dz)
         ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !> if driftkinetic_implicit=T, then only want to treat vpar . grad (<phi>-phi)*F0 term explicitly;
         !> in this case, zero out dg/dz term (or d(g/F)/dz for full-flux-surface)
         if (driftkinetic_implicit) then
            g0 = 0.
         else
            !> compute dg/dz in k-space and store in g0
            call get_dgdz_variable(g(:, :, :, :, ivmu), ivmu, g0)
            !> if simulating a full flux surface, need to obtain the contribution from parallel streaming
            !> in y-space, so FFT d(g/F)/dz from ky to y
            if (full_flux_surface) then
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                     call swap_kxky(g0(:, :, iz, it), g0_swap)
                     call transform_ky2y(g0_swap, g0y(:, :, iz, it))
                  end do
               end do
            end if
            ! ! if simulating a full flux surface, must calculate F * d/dz (g/F) rather than dg/dz
            ! ! since F=F(y) in this case, avoid multiple Fourier transforms by applying chain rule
            ! ! to z derivative: F * d/dz (g/F) = dg/dz - g * d ln F / dz = dg/dz + g * mu/T * dB/dz
            ! if (full_flux_surface) then
            !    ! transform g and dg/dz from ky to y space and store in g0y and g1y, respectively
            !    g1y = g(:,:,:,:,ivmu)
            !    do it = 1, ntubes
            !       do iz = -nzgrid, nzgrid
            !          call transform_ky2y (g1y(:,:,iz,it), g0y(:,:,iz,it))
            !          ! no longer need g1y so re-use as FFT of dg/dz (g0)
            !          call transform_ky2y (g0(:,:,iz,it), g1y(:,:,iz,it))
            !       end do
            !    end do
            !    ! overwrite g0y with dg/dz + g * mu/T * dB/dz
            !    g0y = g1y + 2.0*mu(imu)*spread(spread(dBdzed,2,nakx),4,ntubes) * g0y
            !    ! g1y no longer needed so can over-write with d<phi>/dz below
            ! end if
         end if

         if (full_flux_surface) then
            !> transform d<chi>/dz (fully explicit) or d(<phi>-phi)/dz (if driftkinetic_implicit)
            !> from kalpha (ky) to alpha (y) space and store in g1y
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call swap_kxky(dgchi_dz(:, :, iz, it), g0_swap)
                  call transform_ky2y(g0_swap, g1y(:, :, iz, it))
               end do
            end do
            !> over-write g0y with d/dz (g/F) + Ze/T * d<phi>/dz (or <phi>-phi for driftkinetic_implicit).
            g0y(:, :, :, :) = g0y(:, :, :, :) + g1y(:, :, :, :) * spec(is)%zt
            ! ! over-write g0y with F * d/dz (g/F) + ZeF/T * d<phi>/dz (or <phi>-phi for driftkinetic_implicit).
            ! g0y(:,:,:,:) = g0y(:,:,:,:) + g1y(:,:,:,:)*spec(is)%zt*maxwell_fac(is) &
            !      * maxwell_vpa(iv,is)*spread(spread(maxwell_mu(:,:,imu,is),2,nakx),4,ntubes)*maxwell_fac(is)

            !> multiply d(g/F)/dz and d<phi>/dz terms with vpa*(b . grad z) and add to source (RHS of GK equation)
            call add_stream_term_ffs(g0y, ivmu, gout(:, :, :, :, ivmu))
         else
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! !!! Bob: commented out for h source term
            ! ia = 1
            ! g0(:, :, :, :) = g0(:, :, :, :) + dgchi_dz(:, :, :, :) * spec(is)%zt * maxwell_fac(is) &
            !                  * maxwell_vpa(iv, is) * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx), 4, ntubes)
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! multiply dg/dz with vpa*(b . grad z) and add to source (RHS of GK equation)
            call add_stream_term(g0, ivmu, gout(:, :, :, :, ivmu))
         end if

      end do

      !> deallocate intermediate arrays used in this subroutine
      deallocate (g0, dgchi_dz)
      if (full_flux_surface) deallocate (g0y, g1y, g0_swap)

      !> finish timing the subroutine
      if (proc0) call time_message(.false., time_parallel_streaming(:, 1), ' Stream advance')

   end subroutine advance_parallel_streaming_explicit

   subroutine add_parallel_streaming_radial_variation(g, gout, rhs)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use job_manage, only: time_message
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use species, only: spec
      use gyro_averages, only: gyro_average, gyro_average_j1
      use fields_arrays, only: phi, phi_corr_QN, phi_corr_GA

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      !the next input/output is for quasineutrality and gyroaveraging corrections
      !that go directly in the RHS (since they don't require further FFTs)
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: rhs

      integer :: ivmu, iv, imu, is, it, ia, iz

      complex, dimension(:, :, :, :), allocatable :: g0, g1, g2, g3
      complex, dimension(:, :), allocatable :: g0k

      allocate (g0(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g1(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g2(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g3(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (g0k(naky, nakx)); g0k = 0

      ! parallel streaming stays in ky,kx,z space with ky,kx,z local

      ia = 1
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! obtain <phi>
         ! get d<phi>/dz, with z the parallel coordinate and store in g1
         call gyro_average(phi, ivmu, g0)
         call get_dgdz_centered(g0, ivmu, g1)

         ! get variation in gyroaveraging and store in g2
         call get_dgdz_centered(phi_corr_GA(:, :, :, :, ivmu), ivmu, g2)

         ! get variation in quasineutrality and store in g3
         call gyro_average(phi_corr_QN, ivmu, g0)
         call get_dgdz_centered(g0, ivmu, g3)

         call get_dgdz(g(:, :, :, :, ivmu), ivmu, g0)

         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid

    !!#1 - variation in gradpar
               g0k = g0(:, :, iz, it) &
                     + g1(:, :, iz, it) * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)

               g0k = g0k * stream_rad_var1(iz, iv, is)

    !!#2 - variation in F_s/T_s
               g0k = g0k + g1(:, :, iz, it) * stream_rad_var2(ia, iz, ivmu)

               gout(:, :, iz, it, ivmu) = gout(:, :, iz, it, ivmu) + g0k

    !!#3 - variation in the gyroaveraging and quasineutrality of phi
    !!     These variations already have the linear part calculated, so
    !!     ad it into the rhs directly
               g0k = spec(is)%zt * stream(1, iz, iv, is) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is) &
                     * (g2(:, :, iz, it) + g3(:, :, iz, it))

               rhs(:, :, iz, it, ivmu) = rhs(:, :, iz, it, ivmu) + g0k

            end do
         end do
      end do
      deallocate (g0, g1, g2, g3, g0k)

   end subroutine add_parallel_streaming_radial_variation

   subroutine get_dgdz_variable(g, ivmu, dgdz)

      use run_parameters, only: zed_upwind_explicit
      use kt_grids, only: nakx, naky
      use zgrid, only: nzgrid, ntubes

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, it, iky, iv
      complex, dimension(2) :: gleft, gright
      complex, dimension(:, :, :, :), allocatable :: dgdz_tmp

      allocate (dgdz_tmp(naky, nakx, -nzgrid:nzgrid, ntubes))
      call get_dgdz(g, ivmu, dgdz_tmp)        ! Third order upwind derivative
      call get_dgdz_centered(g, ivmu, dgdz)   ! Second order centered derivative

      dgdz = (1 - zed_upwind_explicit) * dgdz + zed_upwind_explicit * dgdz_tmp

      deallocate (dgdz_tmp)
   end subroutine get_dgdz_variable

   subroutine get_dgdz(g, ivmu, dgdz)

      use finite_differences, only: third_order_upwind_zed
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx
      use zgrid, only: nzgrid, delzed, ntubes
      use extended_zgrid, only: neigen, nsegments
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: fill_zed_ghost_zones
      use extended_zgrid, only: periodic
      use kt_grids, only: naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, it, iky, iv
      complex, dimension(2) :: gleft, gright

      ! FLAG -- assuming delta zed is equally spaced below!
      iv = iv_idx(vmu_lo, ivmu)
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :), gleft, gright)
                  ! now get dg/dz
                  call third_order_upwind_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                                              g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                              delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                                              dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dgdz

   subroutine get_dgdz_centered(g, ivmu, dgdz)

      use finite_differences, only: second_order_centered_zed
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx
      use zgrid, only: nzgrid, delzed, ntubes
      use extended_zgrid, only: neigen, nsegments
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: fill_zed_ghost_zones
      use extended_zgrid, only: periodic
      use kt_grids, only: naky

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz
      integer, intent(in) :: ivmu

      integer :: iseg, ie, iky, iv, it
      complex, dimension(2) :: gleft, gright
      ! FLAG -- assuming delta zed is equally spaced below!
      iv = iv_idx(vmu_lo, ivmu)
      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g(:, :, :, :), gleft, gright)
                  ! now get dg/dz
                  call second_order_centered_zed(iz_low(iseg), iseg, nsegments(ie, iky), &
                                                 g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                                 delzed(0), stream_sign(iv), gleft, gright, periodic(iky), &
                                                 dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do
   end subroutine get_dgdz_centered

! subroutine get_dgdz_variable (g, ivmu, dgdz)

!    use finite_differences, only: fd_variable_upwinding_zed
!    use stella_layouts, only: vmu_lo
!    use stella_layouts, only: iv_idx
!    use zgrid, only: nzgrid, delzed, ntubes
!    use extended_zgrid, only: neigen, nsegments
!    use extended_zgrid, only: iz_low, iz_up
!    use extended_zgrid, only: ikxmod
!    use extended_zgrid, only: fill_zed_ghost_zones
!    use extended_zgrid, only: periodic
!    use run_parameters, only: zed_upwind
!    use kt_grids, only: naky

!    implicit none

!    complex, dimension (:,:,-nzgrid:,:), intent (in) :: g
!    complex, dimension (:,:,-nzgrid:,:), intent (out) :: dgdz
!    integer, intent (in) :: ivmu

!    integer :: iseg, ie, iky, iv, it
!    complex, dimension (2) :: gleft, gright
!    ! FLAG -- assuming delta zed is equally spaced below!
!     iv = iv_idx(vmu_lo,ivmu)
!     do iky = 1, naky
!       do it = 1, ntubes
!         do ie = 1, neigen(iky)
!           do iseg = 1, nsegments(ie,iky)
!              ! first fill in ghost zones at boundaries in g(z)
!              call fill_zed_ghost_zones (it, iseg, ie, iky, g(:,:,:,:), gleft, gright)
   ! now get dg/dz
!              call fd_variable_upwinding_zed (iz_low(iseg), iseg, nsegments(ie,iky), &
!                   g(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it), &
!                   delzed(0), stream_sign(iv), zed_upwind,gleft, gright, periodic(iky), &
!                   dgdz(iky,ikxmod(iseg,ie,iky),iz_low(iseg):iz_up(iseg),it))
!           end do
!         end do
!       enddo
!     end do
! end subroutine get_dgdz_variable

   subroutine add_stream_term(g, ivmu, src)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use zgrid, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src
      integer, intent(in) :: ivmu

      integer :: iz, iv, is

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      do iz = -nzgrid, nzgrid
         src(:, :, iz, :) = src(:, :, iz, :) + stream(1, iz, iv, is) * g(:, :, iz, :)
      end do

   end subroutine add_stream_term

   subroutine add_stream_term_ffs(g, ivmu, src)

      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use zgrid, only: nzgrid
      use kt_grids, only: ny

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: src
      integer, intent(in) :: ivmu

      integer :: iz, iy, iv, is

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      do iz = -nzgrid, nzgrid
         do iy = 1, ny
            src(iy, :, iz, :) = src(iy, :, iz, :) + stream(iy, iz, iv, is) * g(iy, :, iz, :)
         end do
      end do

   end subroutine add_stream_term_ffs

   subroutine advance_parallel_streaming_implicit(g, phi, apar, bpar)

      use mp, only: proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use dist_fn_arrays, only: g1, h
      use run_parameters, only: stream_matrix_inversion
      use fields, only: advance_fields, fields_updated
      use fields, only: get_h, get_gbar

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar

      integer :: ivmu
      complex, dimension(:, :, :, :), allocatable :: deltaphi, deltaapar, deltabpar

      if (proc0) call time_message(.false., time_parallel_streaming(:, 1), ' Stream advance')
      fields_updated = .false.

      allocate (deltaphi(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (deltaapar(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (deltabpar(naky, nakx, -nzgrid:nzgrid, ntubes))

      call advance_fields(g, phi, apar, bpar, dist='gbar')

      ! Get h^{n}
      call get_h(g, phi, apar, bpar, h)

      ! save the incoming h, as they will be needed later
      ! Store in the variable g1, for historical reasons
      g1 = h

      if (proc0) call time_message(.false., time_parallel_streaming(:, 2), ' (bidiagonal solve)')
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! obtain RHS of inhomogeneous GK eqn;
         ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g_{inh}^{n+1}
         ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n}
         ! + (1-alph)/2*dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d<phi^{n}>/dz
         !
         ! New inhomogeneous equation looks like:
         ! (1+(1+alph)/2*dt*vpa*gradpar*d/dz)h_{inh}^{n+1}
         ! = g^{n} - dt*vpa*gradpar*((1-alph)/2)*dh^{n}/dz
         call get_gke_rhs(ivmu, g1(:, :, :, :, ivmu), &
                          phi, apar, bpar, deltaphi, deltaapar, deltabpar, &
                          h(:, :, :, :, ivmu), eqn='inhomogeneous')

         if (stream_matrix_inversion) then
            ! solve (I + (1+alph)/2*dt*vpa . grad)h_{inh}^{n+1} = RHS
            ! h = RHS is input and overwritten by h = h_{inh}^{n+1}
            call invert_parstream(ivmu, h(:, :, :, :, ivmu))
         else
            call sweep_g_zed(ivmu, h(:, :, :, :, ivmu))
         end if
      end do
      if (proc0) call time_message(.false., time_parallel_streaming(:, 2), ' (bidiagonal solve)')

      fields_updated = .false.

      ! we now have h_{inh}^{n+1}
      ! calculate associated fields (phi_{inh}^{n+1}, apar_{inh}^{n+1}, bpar_{inh}^{n+1})
      ! and store in deltaphi etc.
      call advance_fields(h, deltaphi, deltaapar, deltabpar, dist='h')

      ! Now calculate the inhomogeneous change in the fields, delta_fields_{inh}
      deltaphi = deltaphi - phi
      deltaapar = deltaapar - apar
      deltabpar = deltabpar - bpar

      ! solve response_matrix*fields^{n+1} = delta_fields_{inh}
      ! where fields=(phi,apar,bpar)
      ! deltaphi, deltaapar, deltabpar = delta_phi_{inh}, delta_apar_{inh}, delta_bpar_{inh}
      ! is input and overwritten by
      ! deltaphi, deltaapar, deltabpar = delta_phi, delta_apar, delta_bpar
      !                                = phi^{n+1} - phi^{n}, ...
      if (proc0) call time_message(.false., time_parallel_streaming(:, 3), ' (back substitution)')
      call invert_parstream_response(deltaphi, deltaapar, deltabpar)
      if (proc0) call time_message(.false., time_parallel_streaming(:, 3), ' (back substitution)')

      if (proc0) call time_message(.false., time_parallel_streaming(:, 2), ' (bidiagonal solve)')
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         ! now have phi^{n+1} for non-negative kx
         ! obtain RHS of GK eqn;
         ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1}
         ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n}
         ! + dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d/dz((1+alph)/2*<phi^{n+1}>+(1-alph)/2*<phi^{n}>)
         !
         ! New equation looks like:
         ! (1+(1+alph)/2*dt*vpa*gradpar*d/dz)h^{n+1}
         ! = g^{n} + Z/T exp(-v^2) <chi^{n+1}> - dt*vpa*gradpar*((1-alph)/2)*dh^{n}/dz
         call get_gke_rhs(ivmu, g1(:, :, :, :, ivmu), &
                          phi, apar, bpar, deltaphi, deltaapar, deltabpar, &
                          h(:, :, :, :, ivmu), eqn='full')

         if (stream_matrix_inversion) then
            ! solve (1+(1+alph)/2*dt*vpa*gradpar*d/dz)h^{n+1} = RHS
            ! h = RHS is input and overwritten by h = h^{n+1}
            call invert_parstream(ivmu, h(:, :, :, :, ivmu))
         else
            call sweep_g_zed(ivmu, h(:, :, :, :, ivmu))
         end if
      end do
      if (proc0) call time_message(.false., time_parallel_streaming(:, 2), ' (bidiagonal solve)')

      fields_updated = .false.
      call advance_fields(h, phi, apar, bpar, dist="h")
      ! Calculate g^{n+1} = h^{n+1} + <chi^{n+1}>
      call get_gbar(h, phi, apar, bpar, g)
      if (proc0) call time_message(.false., time_parallel_streaming(:, 1), ' Stream advance')

      deallocate (deltaphi)
      deallocate (deltaapar)
      deallocate (deltabpar)

   end subroutine advance_parallel_streaming_implicit

   !> Get the RHS of the parallel streaming piece of the GKE,
   !> when using implicit scheme. Optionally, this can also include the
   !> wdrift and wstar terms (if these are treated implicitly.)
   !> h is output, corresponding to RHS
   !> RHS = RHS_stream + RHS_drifts
   !> RHS_stream = h^{n} + Z/T exp(-v^2) <chi^{n+1}-chi^{n}> - dt*vpa*gradpar*((1-alph)/2)*dh^{n}/dz
   !> RHS_drifts = (1-alph)/2)*(dh^{n}/dx*wdriftx + dh^{n}/dy*wdrifty)
   !>              + wstar*dchi^{n}/dy + (1+alph)/2)*wstar* d/dy(chi^{n+1}-chi^{n})
   subroutine get_gke_rhs(ivmu, hold, phiold, aparold, bparold, &
                          deltaphi, deltaapar, deltabpar, h, eqn)

      use mp, only: mp_abort
      use stella_time, only: code_dt
      use zgrid, only: nzgrid, ntubes
      use species, only: spec
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use kt_grids, only: naky, nakx
      use gyro_averages, only: gyro_average
      use fields, only: get_gyroaverage_chi, get_chi, get_dchidy
      use vpamu_grids, only: vpa, mu
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use stella_geometry, only: dbdzed
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dfneo_dvpa
      use run_parameters, only: time_upwind
      use run_parameters, only: driftkinetic_implicit
      use run_parameters, only: maxwellian_inside_zed_derivative
      use run_parameters, only: stream_drifts_implicit
      use spectral_derivatives, only: get_dgdy, get_dgdx
      use dist_fn_arrays, only: wstar, wdrifty_g, wdriftx_g

      implicit none

      integer, intent(in) :: ivmu
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: hold
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: phiold, deltaphi
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: aparold, deltaapar
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: bparold, deltabpar
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: h
      character(*), intent(in) :: eqn

      integer :: iv, imu, is, iz, ia
      real :: tupwnd1, tupwnd2, inhomogeneous_fac, homogeneous_fac, fac
      ! real, dimension(:), allocatable :: vpadf0dE_fac
      real, dimension(:), allocatable :: gp
      real, dimension(:), allocatable :: maxwell_mu_centered
      real, dimension(:), allocatable :: wdrifty_g_centered, wdriftx_g_centered, wstar_centered
      complex, dimension(:, :, :, :), allocatable :: dhdz
      complex, dimension(:, :, :, :), allocatable :: deltachi
      complex, dimension(:, :, :, :), allocatable :: dhdy, dhdx, dchiolddy, ddeltachidy

      ! allocate (vpadf0dE_fac(-nzgrid:nzgrid))
      allocate (maxwell_mu_centered(-nzgrid:nzgrid))
      allocate (gp(-nzgrid:nzgrid))
      allocate (dhdz(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (deltachi(naky, nakx, -nzgrid:nzgrid, ntubes))

      ia = 1

      tupwnd1 = 0.5 * (1.0 - time_upwind)
      tupwnd2 = 0.5 * (1.0 + time_upwind)

      ! Set flags so that certain terms are ignored when the equation isn't
      ! "full". (The inhomogeneous and full equations are solved every timestep,
      ! the homogeneous equation only computed when the response matrices are being
      ! constructed.
      ! Could refactor to avoid calculating source terms unnecessarily.
      if (eqn == 'full') then
         homogeneous_fac = 1.0
         inhomogeneous_fac = 1.0
      else if (eqn == 'inhomogeneous') then
         homogeneous_fac = 0.0
         inhomogeneous_fac = 1.0
      else if (eqn == 'homogeneous') then
         homogeneous_fac = 1.0
         inhomogeneous_fac = 0.0
      else
         call mp_abort("Flavour of streaming equation not recognised")
      end if

      !! Old comments - delete?
      ! now have phi^{n+1} for non-negative kx
      ! obtain RHS of GK eqn;
      ! i.e., (1+(1+alph)/2*dt*vpa*gradpar*d/dz)g^{n+1}
      ! = (1-(1-alph)/2*dt*vpa*gradpar*d/dz)g^{n}
      ! + dt*Ze*dlnF0/dE*exp(-vpa^2)*vpa*b.gradz*d/dz((1+alph)/2*<phi^{n+1}>+(1-alph)/2*<phi^{n}>
      ! + drift source term

      ! Construct the RHS; the equation looks like:
      ! LHS = (1+(1+alph)/2*dt*vpa*gradpar*d/dz)h^{n+1}
      ! RHS = h^{n} - Z/T exp(-v^2) <delta chi^{n+1}> - dt*vpa*gradpar*((1-alph)/2)*dh^{n}/dz
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! obtain dh^{n}/dz and store in dhdz
      ! NB: could eliminate this calculation at the expense of memory
      ! as this was calculated previously
      call get_dzed(iv, hold, dhdz)

      ! set g to be chi or <chi> depending on whether parallel streaming is
      ! implicit or only implicit in the kperp = 0 (drift kinetic) piece
      if (driftkinetic_implicit) then
         call get_chi(ivmu, deltaphi, deltaapar, deltabpar, deltachi)
         ! call get_chi(ivmu, phiold, aparold, bpar, chiold)
      else
         call get_gyroaverage_chi(ivmu, deltaphi, deltaapar, deltabpar, deltachi)
      end if

      if (maxwellian_inside_zed_derivative) then
         call mp_abort("Maxwellian inside zed derivative currently not supported")
         !    ! From gbar source
         !    ! ! obtain d(exp(-mu*B/T)*<phi>)/dz and store in dchidz
         !    ! g = g * spread(spread(spread(maxwell_mu(ia, :, imu, is), 1, naky), 2, nakx), 4, ntubes)
         !    ! call get_dzed(iv, g, dchidz)
         !    ! ! get <phi>*exp(-mu*B/T)*dB/dz at cell centres
         !    ! g = g * spread(spread(spread(dbdzed(ia, :), 1, naky), 2, nakx), 4, ntubes)
         !    ! call center_zed(iv, g)
         !    ! ! construct d(<phi>*exp(-mu*B/T))/dz + 2*mu*<phi>*exp(-mu*B/T)*dB/dz
         !    ! ! = d<phi>/dz * exp(-mu*B/T)
         !    ! dchidz = dchidz + 2.0 * mu(imu) * g
      else
         ! center Maxwellian factor in mu
         maxwell_mu_centered = maxwell_mu(ia, :, imu, is)
         call center_zed_midpoint(iv, maxwell_mu_centered)
      end if

      ! ! NB: could do this once at beginning of simulation to speed things up
      ! ! this is vpa*Z/T*exp(-vpa^2)
      ! vpadf0dE_fac = vpa(iv) * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_fac(is)
      ! ! if including neoclassical correction to equilibrium distribution function
      ! ! then must also account for -vpa*dF_neo/dvpa*Z/T
      ! ! CHECK TO ENSURE THAT DFNEO_DVPA EXCLUDES EXP(-MU*B/T) FACTOR !!
      if (include_neoclassical_terms) then
         call mp_abort("include_neoclassical_terms currently not supported")
         ! do iz = -nzgrid, nzgrid
         !    vpadf0dE_fac(iz) = vpadf0dE_fac(iz) - 0.5 * dfneo_dvpa(1, iz, ivmu) * spec(is)%zt
         ! end do
         ! call center_zed(iv, vpadf0dE_fac)
      end if

      ! Center h and delta chi in zed.
      h = hold
      call center_zed(iv, h)
      call center_zed(iv, deltachi)

      if (stream_sign(iv) > 0) then
         gp = gradpar_c(:, -1)
      else
         gp = gradpar_c(:, 1)
      end if

      if (stream_drifts_implicit) then
         ! The drift source term is
         ! RHS = drifts_inh + drifts_hom
         ! drifts_inh = + (1+u_t)/2 * ( (wdrift_x * dchi^n/dy) + (wdrift_y * dchi^n/dx) )
         !              + (wstar * dchi^n/dy)
         ! drifts_hom = + (1+u_t)/2 * wstar * d/dy(delta chi)
         !
         ! Center the prefactors. Could save computational cost at the expense of
         ! memory by doing this once at the start.

         ! Allocate memory
         allocate (wdrifty_g_centered(-nzgrid:nzgrid))
         allocate (wdriftx_g_centered(-nzgrid:nzgrid))
         allocate (wstar_centered(-nzgrid:nzgrid))
         allocate (dhdy(naky, nakx, -nzgrid:nzgrid, ntubes))
         allocate (dhdx(naky, nakx, -nzgrid:nzgrid, ntubes))
         allocate (dchiolddy(naky, nakx, -nzgrid:nzgrid, ntubes))
         allocate (ddeltachidy(naky, nakx, -nzgrid:nzgrid, ntubes))

         ! Center the variables in z.
         wdrifty_g_centered = wdrifty_g(1, :, ivmu)
         call center_zed(iv, wdrifty_g_centered)
         wdriftx_g_centered = wdriftx_g(1, :, ivmu)
         call center_zed(iv, wdriftx_g_centered)
         wstar_centered = wstar(1, :, ivmu)
         call center_zed(iv, wstar_centered)

         call get_dgdy(hold, dhdy)
         call center_zed(iv, dhdy)
         call get_dgdx(hold, dhdx)
         call center_zed(iv, dhdx)
         call get_dchidy(ivmu, phiold, aparold, bparold, dchiolddy)
         call center_zed(iv, dchiolddy)
         call get_dchidy(ivmu, deltaphi, deltaapar, deltabpar, ddeltachidy)
         call center_zed(iv, ddeltachidy)

      end if
      ! construct RHS of GK eqn
      ! RHS = h^{n} - Z/T exp(-v^2) <delta chi^{n+1}> - dt*vpa*gradpar*((1-alph)/2)*dh^{n}/dz
      ! where h^{n} and <delta chi> are centered in zed
      fac = code_dt * spec(is)%stm_psi0
      do iz = -nzgrid, nzgrid
         h(:, :, iz, :) = inhomogeneous_fac * h(:, :, iz, :) &
                          + homogeneous_fac * spec(is)%zt * maxwell_vpa(iv, is) * maxwell_mu_centered(iz) * maxwell_fac(is) * deltachi(:, :, iz, :) &
                          - inhomogeneous_fac * fac * gp(iz) * tupwnd1 * vpa(iv) * dhdz(:, :, iz, :)

         if (stream_drifts_implicit) then
            h(:, :, iz, :) = h(:, :, iz, :) &
                             + inhomogeneous_fac * tupwnd1 * dhdy(:, :, iz, :) * wdrifty_g_centered(iz) & ! inh. wdrifty term
                             + inhomogeneous_fac * tupwnd1 * dhdx(:, :, iz, :) * wdriftx_g_centered(iz) & ! inh. wdriftx term
                             + inhomogeneous_fac * dchiolddy(:, :, iz, :) * wstar_centered(iz) & ! inh. wstar term
                             + homogeneous_fac * tupwnd2 * ddeltachidy(:, :, iz, :) * wstar_centered(iz)  ! hom. wstar term
         end if
      end do

      deallocate (maxwell_mu_centered)
      deallocate (gp)
      deallocate (dhdz)
      deallocate (deltachi)
      if (allocated(wdrifty_g_centered)) deallocate (wdrifty_g_centered)
      if (allocated(wdriftx_g_centered)) deallocate (wdriftx_g_centered)
      if (allocated(wstar_centered)) deallocate (wstar_centered)
      if (allocated(dhdy)) deallocate (dhdy)
      if (allocated(dhdx)) deallocate (dhdx)
      if (allocated(dchiolddy)) deallocate (dchiolddy)
      if (allocated(ddeltachidy)) deallocate (ddeltachidy)

   end subroutine get_gke_rhs

   ! solve (I + (1+alph)/2*dt*vpa . grad)g^{n+1} = RHS
   ! g = RHS is input and overwritten by g = g^{n+1}
   subroutine invert_parstream(ivmu, g)

      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: neigen
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use extended_zgrid, only: map_to_extended_zgrid
      use extended_zgrid, only: map_from_extended_zgrid
      use extended_zgrid, only: periodic
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use kt_grids, only: naky

      implicit none

      integer, intent(in) :: ivmu
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: g

      integer :: iv, is
      integer :: iky, ie, it
      integer :: ulim, sgn
      complex, dimension(:), allocatable :: gext

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      sgn = stream_sign(iv)

      do iky = 1, naky
         if (periodic(iky)) then
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  call sweep_zed_zonal(iky, iv, is, sgn, g(iky, ie, :, it))
               end do
            end do
         else
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  allocate (gext(nsegments(ie, iky) * nzed_segment + 1))
                  ! get g on extended domain in zed
                  call map_to_extended_zgrid(it, ie, iky, g(iky, :, :, :), gext, ulim)
                  ! solve (I + (1+alph)/2*dt*vpa . grad)g_{inh}^{n+1} = RHS
                  call stream_tridiagonal_solve(iky, ie, ivmu, gext(:ulim))
                  ! extract g from extended domain in zed
                  call map_from_extended_zgrid(it, ie, iky, gext, g(iky, :, :, :))
                  deallocate (gext)
               end do
            end do
         end if
      end do

   end subroutine invert_parstream

   subroutine stream_tridiagonal_solve(iky, ie, ivmu, g)

      use constants, only: zi
      use finite_differences, only: tridag
      use extended_zgrid, only: iz_low, iz_up, ikxmod
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use dist_fn_arrays, only: wdrifty_g, wdriftx_g
      use run_parameters, only: stream_drifts_implicit
      use stella_layouts, only: iv_idx, is_idx, vmu_lo
      use kt_grids, only: akx, aky

      implicit none

      integer, intent(in) :: iky, ie, ivmu
      complex, dimension(:), intent(in out) :: g

      integer :: iseg, llim, ulim, n
      integer :: nz, nseg_max, sgn, n_ext
      integer :: ia, ikx
      integer :: iv, is
      complex, dimension(:), allocatable :: a, b, c

      ia = 1

      ! avoid double-counting at boundaries between 2pi segments
      nz = nzed_segment
      nseg_max = nsegments(ie, iky)

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      sgn = stream_sign(iv)

      n_ext = nseg_max * nz + 1
      allocate (a(n_ext))
      allocate (b(n_ext))
      allocate (c(n_ext))

      iseg = 1
      llim = 1; ulim = nz + 1
      a(llim:ulim) = stream_tri_a1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_a2(llim:ulim, sgn)
      b(llim:ulim) = stream_tri_b1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_b2(llim:ulim, sgn)
      c(llim:ulim) = stream_tri_c1(llim:ulim, sgn) &
                     - stream(ia, iz_low(iseg):iz_up(iseg), iv, is) * stream_tri_c2(llim:ulim, sgn)

      if (stream_drifts_implicit) then
         ikx = ikxmod(iseg, ie, iky)

         a(llim:ulim) = a(llim:ulim) - drift_tri_a1(llim:ulim, sgn) &
                        * (zi * wdrifty_g(ia, iz_low(iseg):iz_up(iseg), ivmu) * aky(iky) &
                           + zi * wdriftx_g(ia, iz_low(iseg):iz_up(iseg), ivmu) * akx(ikx))

         b(llim:ulim) = b(llim:ulim) - drift_tri_b1(llim:ulim, sgn) &
                        * (zi * wdrifty_g(ia, iz_low(iseg):iz_up(iseg), ivmu) * aky(iky) &
                           + zi * wdriftx_g(ia, iz_low(iseg):iz_up(iseg), ivmu) * akx(ikx))

         c(llim:ulim) = c(llim:ulim) - drift_tri_c1(llim:ulim, sgn) &
                        * (zi * wdrifty_g(ia, iz_low(iseg):iz_up(iseg), ivmu) * aky(iky) &
                           + zi * wdriftx_g(ia, iz_low(iseg):iz_up(iseg), ivmu) * akx(ikx))
      end if

      if (nsegments(ie, iky) > 1) then
         do iseg = 2, nsegments(ie, iky)
            llim = ulim + 1
            ulim = llim + nz - 1
            a(llim:ulim) = stream_tri_a1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_a2(llim:ulim, sgn)
            b(llim:ulim) = stream_tri_b1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_b2(llim:ulim, sgn)
            c(llim:ulim) = stream_tri_c1(llim:ulim, sgn) &
                           - stream(ia, iz_low(iseg) + 1:iz_up(iseg), iv, is) * stream_tri_c2(llim:ulim, sgn)

            if (stream_drifts_implicit) then
               ikx = ikxmod(iseg, ie, iky)
               a(llim:ulim) = a(llim:ulim) - drift_tri_a1(llim:ulim, sgn) &
                              * (zi * wdrifty_g(ia, iz_low(iseg) + 1:iz_up(iseg), ivmu) * aky(iky) &
                                 + zi * wdriftx_g(ia, iz_low(iseg) + 1:iz_up(iseg), ivmu) * akx(ikx))
               b(llim:ulim) = b(llim:ulim) - drift_tri_b1(llim:ulim, sgn) &
                              * (zi * wdrifty_g(ia, iz_low(iseg) + 1:iz_up(iseg), ivmu) * aky(iky) &
                                 + zi * wdriftx_g(ia, iz_low(iseg) + 1:iz_up(iseg), ivmu) * akx(ikx))
               c(llim:ulim) = c(llim:ulim) - drift_tri_c1(llim:ulim, sgn) &
                              * (zi * wdrifty_g(ia, iz_low(iseg) + 1:iz_up(iseg), ivmu) * aky(iky) &
                                 + zi * wdriftx_g(ia, iz_low(iseg) + 1:iz_up(iseg), ivmu) * akx(ikx))
            end if

         end do
      end if
      n = size(stream_tri_a1, 1)
      a(ulim) = stream_tri_a1(n, sgn) - stream(ia, iz_up(nsegments(ie, iky)), iv, is) * stream_tri_a2(n, sgn)
      b(ulim) = stream_tri_b1(n, sgn) - stream(ia, iz_up(nsegments(ie, iky)), iv, is) * stream_tri_b2(n, sgn)
      c(ulim) = 0. ! this line should not be necessary, as c(ulim) should not be accessed by tridag
      if (stream_drifts_implicit) then
         a(ulim) = a(ulim) - drift_tri_a1(n, sgn) &
                   * (zi * wdrifty_g(ia, iz_up(nsegments(ie, iky)), ivmu) * aky(iky) &
                      + zi * wdriftx_g(ia, iz_up(nsegments(ie, iky)), ivmu) * akx(ikx))
         b(ulim) = b(ulim) - drift_tri_b1(n, sgn) &
                   * (zi * wdrifty_g(ia, iz_up(nsegments(ie, iky)), ivmu) * aky(iky) &
                      + zi * wdriftx_g(ia, iz_up(nsegments(ie, iky)), ivmu) * akx(ikx))
      end if
      call tridag(1, a(:ulim), b(:ulim), c(:ulim), g)

      deallocate (a, b, c)

   end subroutine stream_tridiagonal_solve

   ! g= RHS of gke is input
   ! g = g^{n+1} is output
   subroutine sweep_g_zed(ivmu, g)

      use zgrid, only: nzgrid, delzed, ntubes
      use extended_zgrid, only: neigen, nsegments, nzed_segment
      use extended_zgrid, only: map_to_extended_zgrid
      use extended_zgrid, only: map_from_extended_zgrid
      use extended_zgrid, only: periodic
      use kt_grids, only: naky
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, is_idx
      use run_parameters, only: zed_upwind, time_upwind

      implicit none

      integer, intent(in) :: ivmu
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: g

      integer :: iv, is
      integer :: iky, ie, it
      integer :: ulim, sgn
      integer :: iz, izext, iz1, iz2
      real :: fac1, fac2
      complex, dimension(:), allocatable :: gext

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)
      sgn = stream_sign(iv)
      ! will sweep to right (positive vpa) or left (negative vpa)
      ! and solve for g on the extended z-grid
      do iky = 1, naky
         if (periodic(iky)) then
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  call sweep_zed_zonal(iky, iv, is, sgn, g(iky, ie, :, it))
               end do
            end do
         else
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  allocate (gext(nsegments(ie, iky) * nzed_segment + 1))
                  ! get g on extended domain in zed
                  call map_to_extended_zgrid(it, ie, iky, g(iky, :, :, :), gext, ulim)
                  if (sgn < 0) then
                     iz1 = 1; iz2 = ulim
                  else
                     iz1 = ulim; iz2 = 1
                  end if
                  izext = iz1; iz = sgn * nzgrid
                  fac1 = 1.0 + zed_upwind + sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
                  gext(izext) = gext(izext) * 2.0 / fac1
                  do izext = iz1 - sgn, iz2, -sgn
                     if (iz == -sgn * nzgrid) then
                        iz = sgn * nzgrid - sgn
                     else
                        iz = iz - sgn
                     end if
                     fac1 = 1.0 + zed_upwind + sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
                     fac2 = 1.0 - zed_upwind - sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
                     gext(izext) = (-gext(izext + sgn) * fac2 + 2.0 * gext(izext)) / fac1
                  end do
                  ! extract g from extended domain in zed
                  call map_from_extended_zgrid(it, ie, iky, gext, g(iky, :, :, :))
                  deallocate (gext)
               end do
            end do
         end if
      end do

   end subroutine sweep_g_zed

   subroutine sweep_zed_zonal(iky, iv, is, sgn, g)

      use zgrid, only: nzgrid, delzed
      use extended_zgrid, only: phase_shift
      use run_parameters, only: zed_upwind, time_upwind

      implicit none

      integer, intent(in) :: iky, iv, is, sgn
      complex, dimension(-nzgrid:), intent(in out) :: g

      integer :: iz, iz1, iz2
      real :: fac1, fac2
      complex :: pf
      complex, dimension(:), allocatable :: gcf, gpi

      allocate (gpi(-nzgrid:nzgrid))
      allocate (gcf(-nzgrid:nzgrid))
      ! ky=0 is 2pi periodic (no extended zgrid)
      ! decompose into complementary function + particular integral
      ! zero BC for particular integral
      ! unit BC for complementary function (no source)
      if (sgn < 0) then
         iz1 = -nzgrid; iz2 = nzgrid
      else
         iz1 = nzgrid; iz2 = -nzgrid
      end if
      pf = phase_shift(iky)**(-sgn)
      gpi(iz1) = 0.; gcf(iz1) = 1.
      do iz = iz1 - sgn, iz2, -sgn
         fac1 = 1.0 + zed_upwind + sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
         fac2 = 1.0 - zed_upwind - sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
         gpi(iz) = (-gpi(iz + sgn) * fac2 + 2.0 * g(iz)) / fac1
         gcf(iz) = -gcf(iz + sgn) * fac2 / fac1
      end do
      ! g = g_PI + (g_PI(pi)/(1-g_CF(pi))) * g_CF
      g = gpi + (pf * gpi(iz2) / (1.0 - pf * gcf(iz2))) * gcf
      deallocate (gpi, gcf)

   end subroutine sweep_zed_zonal

   subroutine invert_parstream_response(phi, apar, bpar)

      use linear_solve, only: lu_back_substitution
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: neigen
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use extended_zgrid, only: map_to_extended_zgrid
      use extended_zgrid, only: map_from_extended_zgrid
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: periodic, phase_shift
      use run_parameters, only: fphi, fapar, fbpar
      use kt_grids, only: naky
      use fields_arrays, only: response_matrix

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar, bpar

      integer :: iky, ie, it, ulim, nfields, ifield, nresponse_per_field
      integer :: ikx
      complex, dimension(:), allocatable :: fields_ext, phi_ext, apar_ext, bpar_ext

      nfields = 0
      if (fphi > epsilon(0.)) nfields = nfields + 1
      if (fapar > epsilon(0.)) nfields = nfields + 1
      if (fbpar > epsilon(0.)) nfields = nfields + 1

      ! need to put the fields into extended zed grid
      do iky = 1, naky
         ! avoid double counting of periodic endpoints for zonal (and any other periodic) modes
         if (periodic(iky)) then
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  ikx = ikxmod(1, ie, iky)
                  ! One fewer unique point for periodic modes.
                  ! nsegments = 1
                  nresponse_per_field = nsegments(ie, iky) * nzed_segment
                  allocate (fields_ext(nresponse_per_field * nfields))
                  ! put all the fields into a single field: fields_ext=(phi,apar,bpar)
                  ifield = 0
                  if (fphi > epsilon(0.)) then
                     ifield = ifield + 1
                     fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field) = phi(iky, ikx, :nzgrid - 1, it)
                  end if
                  if (fapar > epsilon(0.)) then
                     ifield = ifield + 1
                     fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field) = apar(iky, ikx, :nzgrid - 1, it)
                  end if
                  if (fbpar > epsilon(0.)) then
                     ifield = ifield + 1
                     fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field) = bpar(iky, ikx, :nzgrid - 1, it)
                  end if

                  call lu_back_substitution(response_matrix(iky)%eigen(ie)%zloc, &
                                            response_matrix(iky)%eigen(ie)%idx, fields_ext)
                  ! Get fields on extended grid -> phi, apar, bpar on regular grids
                  ifield = 0
                  if (fphi > epsilon(0.)) then
                     ifield = ifield + 1
                     phi(iky, ikx, :nzgrid - 1, it) = fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field)
                     phi(iky, ikx, nzgrid, it) = phi(iky, ikx, -nzgrid, it) / phase_shift(iky)
                  end if
                  if (fapar > epsilon(0.)) then
                     ifield = ifield + 1
                     apar(iky, ikx, :nzgrid - 1, it) = fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field)
                     apar(iky, ikx, nzgrid, it) = apar(iky, ikx, -nzgrid, it) / phase_shift(iky)
                  end if
                  if (fbpar > epsilon(0.)) then
                     ifield = ifield + 1
                     bpar(iky, ikx, :nzgrid - 1, it) = fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field)
                     bpar(iky, ikx, nzgrid, it) = bpar(iky, ikx, -nzgrid, it) / phase_shift(iky)
                  end if

                  deallocate (fields_ext)
               end do
            end do
         else
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  nresponse_per_field = nsegments(ie, iky) * nzed_segment + 1
                  allocate (fields_ext(nresponse_per_field * nfields))
                  allocate (phi_ext(nresponse_per_field))
                  ifield = 0
                  if (fphi > epsilon(0.)) then
                     ifield = ifield + 1
                     call map_to_extended_zgrid(it, ie, iky, phi(iky, :, :, :), phi_ext, ulim)
                     fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field) = phi_ext
                  end if
                  if (fapar > epsilon(0.)) then
                     ifield = ifield + 1
                     call map_to_extended_zgrid(it, ie, iky, apar(iky, :, :, :), phi_ext, ulim)
                     fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field) = phi_ext
                  end if
                  if (fbpar > epsilon(0.)) then
                     ifield = ifield + 1
                     call map_to_extended_zgrid(it, ie, iky, bpar(iky, :, :, :), phi_ext, ulim)
                     fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field) = phi_ext
                  end if
                  ! solve response_matrix*phi^{n+1} = phi_{inh}^{n+1}
                  call lu_back_substitution(response_matrix(iky)%eigen(ie)%zloc, &
                                            response_matrix(iky)%eigen(ie)%idx, fields_ext)
                  ! Get fields on extended grid -> phi, apar, bpar on regular grids
                  ifield = 0
                  if (fphi > epsilon(0.)) then
                     ifield = ifield + 1
                     phi_ext = fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field)
                     call map_from_extended_zgrid(it, ie, iky, phi_ext, phi(iky, :, :, :))
                  end if
                  if (fapar > epsilon(0.)) then
                     ifield = ifield + 1
                     phi_ext = fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field)
                     call map_from_extended_zgrid(it, ie, iky, phi_ext, apar(iky, :, :, :))
                  end if
                  if (fbpar > epsilon(0.)) then
                     ifield = ifield + 1
                     phi_ext = fields_ext((ifield - 1) * nresponse_per_field + 1:ifield * nresponse_per_field)
                     call map_from_extended_zgrid(it, ie, iky, phi_ext, bpar(iky, :, :, :))
                  end if

                  deallocate (fields_ext)
                  deallocate (phi_ext)
               end do
            end do
         end if
      end do

   end subroutine invert_parstream_response

   subroutine get_dzed(iv, g, dgdz)

      use finite_differences, only: fd_cell_centres_zed
      use kt_grids, only: naky
      use zgrid, only: nzgrid, delzed, ntubes
      use extended_zgrid, only: neigen, nsegments
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: fill_zed_ghost_zones

      implicit none

      integer, intent(in) :: iv
      complex, dimension(:, :, -nzgrid:, :), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(out) :: dgdz

      integer :: iky, ie, iseg, it
      complex, dimension(2) :: gleft, gright

      do it = 1, ntubes
         do iky = 1, naky
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g, gleft, gright)
                  ! get finite difference approximation for dg/dz at cell centres
                  ! iv > nvgrid corresponds to positive vpa, iv <= nvgrid to negative vpa
                  call fd_cell_centres_zed(iz_low(iseg), &
                                           g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                           delzed(0), stream_sign(iv), gleft(2), gright(1), &
                                           dgdz(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
      end do

   end subroutine get_dzed

   subroutine center_zed_extended(iv, g)

      use finite_differences, only: cell_centres_zed
      use kt_grids, only: naky, nakx
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: neigen, nsegments
      use extended_zgrid, only: iz_low, iz_up
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: fill_zed_ghost_zones
      use run_parameters, only: zed_upwind

      implicit none

      integer, intent(in) :: iv
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: g

      integer :: iky, ie, iseg, it
      complex, dimension(2) :: gleft, gright
      complex, dimension(:, :, :), allocatable :: gc

      allocate (gc(nakx, -nzgrid:nzgrid, ntubes))

      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               do iseg = 1, nsegments(ie, iky)
                  ! first fill in ghost zones at boundaries in g(z)
                  call fill_zed_ghost_zones(it, iseg, ie, iky, g, gleft, gright)
                  ! get cell centres values
                  call cell_centres_zed(iz_low(iseg), &
                                        g(iky, ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it), &
                                        zed_upwind, stream_sign(iv), gleft(2), gright(1), &
                                        gc(ikxmod(iseg, ie, iky), iz_low(iseg):iz_up(iseg), it))
               end do
            end do
         end do
         g(iky, :, :, :) = gc
      end do

      deallocate (gc)

   end subroutine center_zed_extended

   subroutine center_zed_segment_real(iv, g)

      use zgrid, only: nzgrid
      use run_parameters, only: zed_upwind

      integer, intent(in) :: iv
      real, dimension(-nzgrid:), intent(in out) :: g

      if (stream_sign(iv) > 0) then
         g(:nzgrid - 1) = 0.5 * ((1.+zed_upwind) * g(:nzgrid - 1) + (1.-zed_upwind) * g(-nzgrid + 1:))
         g(nzgrid) = g(-nzgrid)
      else
         g(-nzgrid + 1:) = 0.5 * ((1.-zed_upwind) * g(:nzgrid - 1) + (1.+zed_upwind) * g(-nzgrid + 1:))
         g(-nzgrid) = g(nzgrid)
      end if

   end subroutine center_zed_segment_real

   subroutine center_zed_midpoint(iv, g)

      use zgrid, only: nzgrid

      integer, intent(in) :: iv
      real, dimension(-nzgrid:), intent(in out) :: g

      if (stream_sign(iv) > 0) then
         g(:nzgrid - 1) = 0.5 * (g(:nzgrid - 1) + g(-nzgrid + 1:))
         g(nzgrid) = g(-nzgrid)
      else
         g(-nzgrid + 1:) = 0.5 * (g(:nzgrid - 1) + g(-nzgrid + 1:))
         g(-nzgrid) = g(nzgrid)
      end if

   end subroutine center_zed_midpoint

   subroutine finish_parallel_streaming

      use run_parameters, only: stream_implicit, driftkinetic_implicit

      implicit none

      if (allocated(stream)) deallocate (stream)
      if (allocated(stream_c)) deallocate (stream_c)
      if (allocated(stream_sign)) deallocate (stream_sign)
      if (allocated(gradpar_c)) deallocate (gradpar_c)
      if (allocated(stream_rad_var1)) deallocate (stream_rad_var1)
      if (allocated(stream_rad_var2)) deallocate (stream_rad_var2)

      if (stream_implicit .or. driftkinetic_implicit) call finish_invert_stream_operator

      parallel_streaming_initialized = .false.

   end subroutine finish_parallel_streaming

   subroutine finish_invert_stream_operator

      implicit none

      if (allocated(stream_tri_a1)) then
         deallocate (stream_tri_a1)
         deallocate (stream_tri_a2)
         deallocate (stream_tri_b1)
         deallocate (stream_tri_b2)
         deallocate (stream_tri_c1)
         deallocate (stream_tri_c2)
      end if
      if (allocated(drift_tri_a1)) then
         deallocate (drift_tri_a1)
         deallocate (drift_tri_b1)
         deallocate (drift_tri_c1)
      end if

   end subroutine finish_invert_stream_operator

end module parallel_streaming
